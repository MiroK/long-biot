# NOTE: we assume the C param in (C + alpha^2/lmbda) to be C = 0 (which
# should be the most difficult case)

from collections import namedtuple
import argparse
from dolfin import *
from petsc4py import PETSc
import sympy as sp
import ulfy
from xii import *

from block.algebraic.petsc import KSP, AMG, LU
from block.block_mat import block_mat
import os

print = PETSc.Sys.Print

BiotParameters = namedtuple('BiotParameters', ('alpha', 'K', 'mu', 'lmbda', 'c'))


from biot3 import (Laplacian, generate_2d_domains, setup_2d_mms, get_system,
                   parse_V_bcs, parse_Q_bcs, get_path)



def get_operator_mixed(boundaries, parameters, *, u_dirichlet_tags, p_dirichlet_tags, W, Wbcs, bdry_tags):
    '''Structure V x (Q x QT) where pressures are coupled'''
    u, p, pT = map(TrialFunction, W)
    v, q, qT = map(TestFunction, W)

    V, Q, QT = W
    
    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                               for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))

    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

    QQQ = FunctionSpace(mesh,
                        MixedElement([Q.ufl_element()]*3))

    p0, p1, p2 = TrialFunctions(QQQ)
    q0, q1, q2 = TestFunctions(QQQ)
                   
    # ---
    a = (c + alpha**2/lmbda)*inner(p0, q0)*dx + inner(K*grad(p0), grad(q0))*dx

    nF, hF = FacetNormal(Q.mesh()), CellDiameter(Q.mesh())
    gammaF = Constant(5)
    # Add Nietsche terms
    for tag in p_dirichlet_tags:
        a += (
            - inner(dot(K*grad(p0), nF), q0)*ds(tag)
            - inner(dot(K*grad(q0), nF), p0)*ds(tag)
            + (K*gammaF/hF)*inner(p0, q0)*ds(tag)
        )
        
    
    a += (-alpha/lmbda)*inner(p1, q0)*dx
    a += (-alpha/lmbda)*inner(p2, q0)*dx

    a += (-alpha/lmbda)*inner(q1, p0)*dx

    a += (1/2/mu + 1/lmbda)*inner(p1, q1)*dx
    a += (1/lmbda)*inner(p2, q1)*dx

    a += (-alpha/lmbda)*inner(q2, p0)*dx    
    a += (1/lmbda)*inner(p1, q2)*dx

    bc_tags = {'T': set(bdry_tags) - set(u_dirichlet_tags),
               'P': set()}
    scale = Constant(1) # FIXME, this will be the thickness
    kappa = alpha**2/(1+lmbda)*scale**2
    k_form, ker = Laplacian((p2, q2), boundaries, bc_tags, kappa=kappa)

    a += k_form + (1/lmbda)*inner(p2, q2)*dx

    # Puttin together
    EE = ii_assemble(a)

    return EE, QQQ

# --------------------------------------------------------------------

if __name__ == '__main__':
    from functools import partial
    import tabulate, argparse    
    import numpy as np
    import os
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Domain
    parser.add_argument('-L', type=int, default=1)
    parser.add_argument('-nrefs', type=int, default=4)
    
    parser.add_argument('-precond', type=str, default='standard')
    
    # bcs
    parser.add_argument('-ubcs', type=str, help='Spec of bcs for momentum: D or T', default='TTDD')
    parser.add_argument('-pbcs', type=str, help='Spec of bcs for mass: P or F', default='FFFF')    
    # Material
    parser.add_argument('-alpha', type=float, default=1E0)
    parser.add_argument('-K', type=float, default=1)
    parser.add_argument('-mu', type=float, default=1)
    parser.add_argument('-lmbda', type=float, default=1)
    parser.add_argument('-c', type=float, default=0)

    args = parser.parse_args()


    result_dir = f'./results/amg_biot3/precond{args.precond}'
    not os.path.exists(result_dir) and os.makedirs(result_dir)

    ignore_keys = ('nrefs', )
    get_path = partial(get_path,
                       args=vars(args), ignore_keys=ignore_keys,
                       result_dir=result_dir)

    mesh_gen = generate_2d_domains(args.L, ncells0=2, kind='crossed')
    boundary_tags = next(mesh_gen)

    #  --------------4------
    # 1|                   |2
    #  --------------3------
    u_dirichlet_tags = parse_V_bcs(args.ubcs, boundary_tags)
    p_dirichlet_tags = parse_Q_bcs(args.pbcs, boundary_tags)


    parameters = BiotParameters(alpha=args.alpha, K=args.K, mu=args.mu, lmbda=args.lmbda, c=args.c)
    
    mms_data = setup_2d_mms(parameters)

    headers = ('h', 'ndofs', 
               'niters', 'lminKSP', 'lmaxKSP', 'condKSP')

    precond = 'hypre'
    length = args.L
    
    h0, errors0, history = None, None, []
    for (level, boundaries) in zip(range(args.nrefs), mesh_gen):
        mesh = boundaries.mesh()
        
        A, b, W, W_bcs = get_system(boundaries, parameters=parameters, data=mms_data,
                                    u_dirichlet_tags=u_dirichlet_tags,
                                    p_dirichlet_tags=p_dirichlet_tags,
                                    bdry_tags=boundary_tags,
                                    pdegrees='2_1_1')

        # [0, 2, 4], [1, 3, 5]
        B, QQ = get_operator_mixed(boundaries, parameters=parameters,
                                   u_dirichlet_tags=u_dirichlet_tags,
                                   p_dirichlet_tags=p_dirichlet_tags,
                                   W=W, Wbcs=W_bcs,
                                   bdry_tags=boundary_tags)

        # [0, 1, 2], [3, 4, 5]
        C, QQ = get_operator(boundaries, parameters=parameters,
                             u_dirichlet_tags=u_dirichlet_tags,
                             p_dirichlet_tags=p_dirichlet_tags,
                             W=W, Wbcs=W_bcs,
                             bdry_tags=boundary_tags)

        # # AA = monolithic(B)
        # #print(AA.array())
        
        # Qi, Qi = QQ
        # elm = MixedElement([Qi.ufl_element()]*2)
        # Z = FunctionSpace(Qi.mesh(), elm)
        # dm = np.array([Z.sub(i).dofmap().dofs() for i in range(2)], dtype='int32')

        
        # perm = PETSc.IS().createGeneral(dm).invertPermutation()
        # D = as_backend_type(monolithic(C)).mat().permute(perm, perm)
        # D.setBlockSize(2)
        # D = PETScMatrix(D)

        # print('>>>', monolithic(monolithic(B) - monolithic(D)).norm('linf'), '<<<')
        
        AA = B
        
        #print(AA.array())
        #from IPython import embed
        #embed()
        # 
        #print((AA-BB).norm('linf'))
        #exit()
        # ---
        solver = PETScKrylovSolver()
        ksp = solver.ksp()

        opts = PETSc.Options()
        opts.setValue('ksp_type', 'cg')
        opts.setValue('ksp_view', None)
        opts.setValue('ksp_monitor_true_residual', None)
        opts.setValue('ksp_rtol', 1E-12)
        opts.setValue('ksp_max_it', 1000)
        # opts.setValue('ksp_norm_type', 'unpreconditioned')
        # opts.setValue('ksp_initial_guess_nonzero', None)
        # opts.setValue('ksp_knoll', None)
        opts.setValue('ksp_converged_reason', None)
        opts.setValue('ksp_view', None)
        opts.setValue('options_view', None)
        
        if precond == 'hypre':
            solver.set_operators(AA, AA)
            opts.setValue('pc_type', 'hypre')
            opts.setValue('pc_hypre_boomeramg_strong_threshold', 0.1)
            opts.setValue('pc_hypre_boomeramg_nodal_coarsen', 1)
            opts.setValue('pc_hypre_boomeramg_vec_interp_variant', 1)
            opts.setValue('pc_hypre_boomeramg_interp_type', 'ext+i')
            opts.setValue('pc_hypre_boomeramg_smooth_type', 'Schwarz-smoothers')

        elif precond == 'gamg':
            solver.set_operators(AA, AA)
            opts.setValue('pc_type', 'gamg')            
            
        elif precond == 'lu':
            pass
            
        else:
            from block.block_base import block_base
            from block.algebraic.petsc.solver import petsc_py_wrapper
            from scipy.sparse import csr_matrix
            import pyamg

            AA_ = as_backend_type(AA).mat()
            indptr, indices, data = AA_.getValuesCSR()
            
            if precond == 'pyamgCSR':
                AA_scipy = csr_matrix((data, indices, indptr), shape=AA_.getSize())
                B = None

            elif precond == 'pyamgBSR':
                AA_scipy = csr_matrix((data, indices, indptr), shape=AA_.getSize()).tobsr((2, 2))
                B = None

            elif precond == 'pyamgBSR':
                AA_scipy = csr_matrix((data, indices, indptr), shape=AA_.getSize()).tobsr((2, 2))

                z = interpolate(Constant((1, 1)), W).vector().get_local()
                z = z/np.linalg.norm(z)
                B = z.reshape((-1, 1))

            else:
                raise ValueError
                
            prec = pyamg.smoothed_aggregation_solver(
                A=AA_scipy,
                B=B,
                aggregate='standard',
                smooth=('energy', {}),
                presmoother=('block_gauss_seidel',
                             {'sweep': 'symmetric'}),
                postsmoother=('block_gauss_seidel',
                              {'sweep': 'symmetric'}),
                improve_candidates=[('block_gauss_seidel',
                                     {'sweep': 'symmetric',
                                      'iterations': 4}),
                                    None],        
            ).aspreconditioner()

            class Precond(block_base):
                def __init__(self, prec):
                    self.prec = prec
                    self.A = AA
                    
                def matvec(self, x):
                    y = x.copy()
                    y.set_local(self.prec@x.get_local())
                    return y
                
                def create_vec(self, which):
                    return PETScVector(AA_.createVecs()[which])

            # ------

            solver.set_operators(AA, AA)
                    
            pc = ksp.getPC()
            pc.setType(PETSc.PC.Type.PYTHON)
            pc.setPythonContext(petsc_py_wrapper(Precond(prec)))

        ksp.setComputeEigenvalues(True)
        ksp.setFromOptions()

        S = StackOperator((2, W[1]))
        b = AA.create_vec()
        b.set_local(np.random.rand(b.size()))

        x = b.copy()
        # bb = monolithic(S*b[1])
        solver.solve(x, b)

        niters = ksp.getIterationNumber()
        
        eigs = ksp.computeEigenvalues()
        lminKSP, lmaxKSP = np.sort(np.abs(eigs))[[0, -1]]
        condKSP = lmaxKSP/lminKSP
        
        # Check convergence
        h = mesh.hmin()
        ndofs = sum(Wi.dim() for Wi in W)

        row = (h, ndofs) + (niters, lminKSP, lmaxKSP, condKSP)
        history.append(row)
        print(tabulate.tabulate(history, headers=headers))

        with open(get_path('cvrg', 'txt'), 'w') as out:
            np.savetxt(out, np.array(history), header=' '.join(headers))
