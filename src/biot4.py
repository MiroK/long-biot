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


from utils import SerializeOperator

print = PETSc.Sys.Print

from biot3 import (BiotParameters, Laplacian, generate_2d_domains, setup_2d_mms,
                   parse_V_bcs, parse_Q_bcs, get_path, SYM)

# ---

def get_system(boundaries, parameters, data, *, u_dirichlet_tags, p_dirichlet_tags, bdry_tags, q_degrees='1_0'):
    '''Three field formulation'''
    mesh = boundaries.mesh()
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    n = FacetNormal(mesh)

    V = VectorFunctionSpace(mesh, 'CG', 2)
    M = FunctionSpace(mesh, 'RT', 1)

    qt_deg, q_deg = map(int, q_degrees.split('_'))
    if qt_deg == 0:
        QT = FunctionSpace(mesh, 'DG', qt_deg)
    else:
        QT = FunctionSpace(mesh, 'CG', qt_deg)
    Q = FunctionSpace(mesh, 'DG', q_deg)
    W = [V, M, QT, Q]
    
    u, l, pT, p = map(TrialFunction, W)
    v, m, qT, q = map(TestFunction, W)

    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                              for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))

    a = block_form(W, 2)
    a[0][0] = inner(2*mu*SYM(grad(u)), SYM(grad(v)))*dx
    a[0][2] = -inner(pT, div(v))*dx
    
    a[1][1] = (1/K)*inner(l, m)*dx
    a[1][3] = -inner(p, div(m))*dx

    a[2][0] = -inner(qT, div(u))*dx
    a[2][2] = (-1/lmbda)*inner(qT, pT)*dx    
    a[2][3] = (alpha/lmbda)*inner(qT, p)*dx

    a[3][1] = -inner(q, div(l))*dx
    a[3][2] = (alpha/lmbda)*inner(pT, q)*dx
    a[3][3] = -(c + alpha**2/lmbda)*inner(p, q)*dx    
    
    u_neumann_tags = bdry_tags - set(u_dirichlet_tags)
    p_neumann_tags = bdry_tags - set(p_dirichlet_tags)

    L = block_form(W, 1)
    
    L[0] = inner(data['f_u'], v)*dx
    L[0] += sum(inner(data['u_neumann'][tag], v)*ds(tag) for tag in u_neumann_tags)
    
    # Pressure eq.
    if p_dirichlet_tags:
        L[1] = -sum(inner(data['p_dirichlet'][tag], dot(m, n))*ds(tag) for tag in p_dirichlet_tags)    

    L[3] = inner(-data['f_p'], q)*dx

    V_bcs = [DirichletBC(V, data['u_dirichlet'][tag], boundaries, tag) for tag in u_dirichlet_tags]
    M_bcs = [DirichletBC(M, data['flux'], boundaries, tag) for tag in p_neumann_tags]
    QT_bcs = []
    Q_bcs = []
    W_bcs = [V_bcs, M_bcs, QT_bcs, Q_bcs]

    A, b = map(ii_assemble, (a, L))
    
    A, b = apply_bc(A, b, bcs=W_bcs)

    return A, b, W, W_bcs


def get_inner_product_standard(boundaries, parameters, *, u_dirichlet_tags, p_dirichlet_tags, W, Wbcs, bdry_tags,
                               inverseQ='lu'):
    '''Structure V x (Q x QT) where pressures are coupled'''
    mesh = boundaries.mesh()
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    n = FacetNormal(mesh)

    V, M, QT, Q = W
    
    u, l, pT, p = map(TrialFunction, W)
    v, m, qT, q = map(TestFunction, W)

    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                              for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))

    a = block_form(W[:2], 2)
    a[0][0] = inner(2*mu*SYM(grad(u)), SYM(grad(v)))*dx
    a[1][1] = (1/K)*inner(l, m)*dx + inner(div(l), div(m))*dx

    BB  = ii_assemble(a)
    BB, _ = apply_bc(BB, b=None, bcs=W_bcs[:2])
    B0, B1 = BB[0][0], BB[1][1]

    QTQ = [QT, Q]

    cA = block_form(QTQ, 2)
    # ---
    cA[0][0] = (1/lmbda + 1/2/mu)*inner(qT, pT)*dx    
    cA[0][1] = (-alpha/lmbda)*inner(qT, p)*dx
    cA[1][0] = (-alpha/lmbda)*inner(pT, q)*dx
    cA[1][1] = (c + alpha**2/lmbda + 1)*inner(p, q)*dx


    cB = block_form(QTQ, 2)
    # ---
    cB[0][0] = (1/lmbda + 1/2/mu)*inner(qT, pT)*dx    
    cB[0][1] = (-alpha/lmbda)*inner(qT, p)*dx
    cB[1][0] = (-alpha/lmbda)*inner(pT, q)*dx
    cB[1][1] = (c + alpha**2/lmbda)*inner(p, q)*dx + inner(K*grad(p), grad(q))*dx

    # Handle DG0 Laplacian
    hF = CentroidDistance(mesh)
    cB[1][1] += K*Constant(1)/avg(hF)*inner(jump(p), jump(q))*dS
    
    # FIXME: add Nitsce bcs?
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    nF, hF = FacetNormal(Q.mesh()), CellDiameter(Q.mesh())
    gammaF = Constant(5)
    # Add Nietsche terms
    for tag in p_dirichlet_tags:
        cB[1][1] += (
            - inner(dot(K*grad(p), nF), q)*ds(tag)
            - inner(dot(K*grad(q), nF), p)*ds(tag)
            + (K*gammaF/hF)*inner(p, q)*ds(tag)
        ) 

    CA, CB = (monolithic(ii_assemble(x)) for x in (cA, cB))
    
    iB0 = LU(B0)
    iB1 = LU(B1)
    iCA, iCB = LU(CA), LU(CB)

    R = ReductionOperator([1, 2, 4], W)
    iBB = R.T*block_diag_mat([iB0, iB1, iCA + iCB])*R

    return BB, iBB


def get_inner_product_espen(boundaries, parameters, *, u_dirichlet_tags, p_dirichlet_tags, W, Wbcs, bdry_tags,
                            inverseQ='lu'):
    '''Structure V x (Q x QT) where pressures are coupled'''
    mesh = boundaries.mesh()
    V, M, QT, Q = W
    
    u, l, pT, p = map(TrialFunction, W)
    v, m, qT, q = map(TestFunction, W)

    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                              for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))

    a = block_form(W[:2], 2)
    a[0][0] = inner(2*mu*SYM(grad(u)), SYM(grad(v)))*dx
    a[1][1] = (1/K)*inner(l, m)*dx + inner(div(l), div(m))*dx

    BB  = ii_assemble(a)
    BB, _ = apply_bc(BB, b=None, bcs=W_bcs[:2])
    B0, B1 = BB[0][0], BB[1][1]

    # -----------------

    if inverseQ == 'amg':
        assert Q.ufl_element() == QT.ufl_element()

    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    nF, hF = FacetNormal(Q.mesh()), CentroidDistance(Q.mesh())
    gammaF = Constant(5)
    # Monolithic operator for stacked pressures    
    QQ = FunctionSpace(mesh,
                       MixedElement([Q.ufl_element(), Q.ufl_element(), QT.ufl_element(), QT.ufl_element()]))
    # ---
    pT, pTb, p, pb = TrialFunctions(QQ)
    qT, qTb, q, qb = TestFunctions(QQ)

    a = (1/2/mu + 1/lmbda)*inner(pT, qT)*dx
    a += (1/lmbda)*inner(pTb, qT)*dx
    a += (-alpha/lmbda)*inner(p, qT)*dx
    a += (-alpha/lmbda)*inner(pb, qT)*dx

    
    a += (1/lmbda)*inner(pT, qTb)*dx
    a += (1/lmbda)*inner(pTb, qTb)*dx
    # FIXME: add Laplacian from Stokes
    bc_tags = {'T': set(bdry_tags) - set(u_dirichlet_tags),
               'P': set()}
    scale = Constant(1) # FIXME, this will be the thickness
    kappa = scale**2/2/mu    
    k_form, ker = Laplacian((pTb, qTb), boundaries, bc_tags, kappa=kappa)
    # ...
    a += k_form    
    
    a += (-alpha/lmbda)*inner(p, qTb)*dx
    a += (-alpha/lmbda)*inner(pb, qTb)*dx

    
    a += (-alpha/lmbda)*inner(pT, q)*dx
    a += (-alpha/lmbda)*inner(pTb, q)*dx
    a += (1 + c + alpha**2/lmbda)*inner(p, q)*dx
    a += (c + alpha**2/lmbda)*inner(pb, q)*dx

    
    a += (-alpha/lmbda)*inner(pT, qb)*dx
    a += (-alpha/lmbda)*inner(pTb, qb)*dx
    a += (c + alpha**2/lmbda)*inner(p, qb)*dx
    a += (c + alpha**2/lmbda)*inner(pb, qb)*dx
    # FIXME: add Laplacian from fluid
    a += inner(K*grad(pb), grad(qb))*dx + K*Constant(1)/avg(hF)*inner(jump(pb), jump(qb))*dS
    # Add Nietsche terms
    for tag in p_dirichlet_tags:
        a += (
            - inner(dot(K*grad(pb), nF), qb)*ds(tag)
            - inner(dot(K*grad(qb), nF), pb)*ds(tag)
            + (K*gammaF/hF)*inner(pb, qb)*ds(tag)
        )

    # -----
    EE = assemble(a)
    
    precond0 = LU(B0)
    precond1 = LU(B1)

    # Putting together
    if inverseQ == 'lu':
        invE = LU(EE)
    elif inverseQ == 'amg':
        invE = AMG(EE,
                   parameters={
                       'pc_hypre_boomeramg_strong_threshold': 0.1,
                       'pc_hypre_boomeramg_nodal_coarsen': 1,
                       'pc_hypre_boomeramg_vec_interp_variant': 1,
                       'pc_hypre_boomeramg_interp_type': 'direct',  #'ext+i',
                       'pc_hypre_boomeramg_smooth_type': 'Schwarz-smoothers',
                       # 'pc_hypre_boomeramg_nodal_relaxation': None,
                   })
    elif inverseQ == 'pyamg':
        from scipy.sparse import csr_matrix
        from block.block_base import block_base
        import pyamg

        EE_ = as_backend_type(EE).mat()
        indptr, indices, data = EE_.getValuesCSR()
        
        EE_scipy = csr_matrix((data, indices, indptr), shape=EE_.getSize()).tobsr((4, 4))

        prec = pyamg.smoothed_aggregation_solver(
            A=EE_scipy,
            aggregate='standard',
            smooth=('energy', {'weighting': 'block'}),
            presmoother=('block_gauss_seidel',
                         {'sweep': 'symmetric'}),
            postsmoother=('block_gauss_seidel',
                          {'sweep': 'symmetric'}),
            improve_candidates=[('block_gauss_seidel',
                                 {'sweep': 'symmetric',
                                  'iterations': 6}),
                                None],        
        ).aspreconditioner()

        class Precond(block_base):
            def __init__(self, prec):
                self.prec = prec
                self.A = EE
                
            def matvec(self, x):
                y = x.copy()
                y.set_local(self.prec@x.get_local())
                return y
            
            def create_vec(self, which):
                return PETScVector(EE_.createVecs()[which])

        invE = Precond(prec)
    # Finally
    precond2 = invE
    
    # And now ...
    S = StackOperator([(2, QT), (2, Q)], pre=[V, M])
    R = ReductionOperator([1, 2, 6], [V, M, QT, QT, Q, Q])

    # Serialization of QQ
    T = SerializeOperator(QQ)

    precond = block_diag_mat([precond0, precond1, T*precond2*T.T])
    
    iBB = S.T*R.T*precond*R*S

    return None, iBB

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
    parser.add_argument('-inverseQ', type=str, default='lu', choices=('lu', 'amg', 'pyamg'))
    # bcs
    parser.add_argument('-ubcs', type=str, help='Spec of bcs for momentum: D or T', default='TTDD')
    parser.add_argument('-pbcs', type=str, help='Spec of bcs for mass: P or F', default='FFFF')    
    # Material
    parser.add_argument('-alpha', type=float, default=1E0)
    parser.add_argument('-K', type=float, default=1)
    parser.add_argument('-mu', type=float, default=0.5)
    parser.add_argument('-lmbda', type=float, default=1)
    parser.add_argument('-c', type=float, default=0.0)

    args = parser.parse_args()


    result_dir = f'./results/biot4/precond{args.precond}'
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


    get_inner_product = {
                         'standard': get_inner_product_standard,
                         'espen': get_inner_product_espen
    }[args.precond]

    parameters = BiotParameters(alpha=args.alpha, K=args.K, mu=args.mu, lmbda=args.lmbda, c=args.c)
    
    mms_data = setup_2d_mms(parameters)

    # Common solver settings
    opts = PETSc.Options()        
    opts.setValue('ksp_rtol', 1E-12)
    opts.setValue('ksp_view_pre', None)
    opts.setValue('ksp_monitor_true_residual', None)
    opts.setValue('ksp_converged_reason', None)
    opts.setValue('options_view', None)

    headers = ('h', 'ndofs', '|GD|', '|GP|',
               '|eu|1', 'reu', '|ez|_div', 'rez', '|epT|0', 'rep', '|ep|0', 'rep0',
               'niters', 'lminKSP', 'lmaxKSP', 'condKSP')

    length = args.L
    
    h0, errors0, history = None, None, []
    for (level, boundaries) in zip(range(args.nrefs), mesh_gen):
        mesh = boundaries.mesh()
        
        A, b, W, W_bcs = get_system(boundaries, parameters=parameters, data=mms_data,
                                    u_dirichlet_tags=u_dirichlet_tags,
                                    p_dirichlet_tags=p_dirichlet_tags,
                                    bdry_tags=boundary_tags,
                                    q_degrees='0_0')

        B, invB = get_inner_product(boundaries, parameters=parameters,
                                    u_dirichlet_tags=u_dirichlet_tags,
                                    p_dirichlet_tags=p_dirichlet_tags,
                                    W=W, Wbcs=W_bcs,
                                    bdry_tags=boundary_tags,
                                    inverseQ=args.inverseQ)

        # Niters
        Ainv = KSP(A, precond=invB, 
                    # PETScOptions
                   ksp_type='minres',
                   ksp_rtol=1E-12,
                   ksp_view=None,
                   ksp_max_it=1_000,
                   ksp_monitor_true_residual=None,
                   ksp_initial_guess_nonzero=1,
                   ksp_converged_reason=None)
        
        x = Ainv*b

        wh = ii_Function(W)
        for i, xi in enumerate(x):
            wh[i].vector().axpy(1, xi)
        niters = len(Ainv.residuals)
        
        eigs = Ainv.eigenvalue_estimates()
        lminKSP, lmaxKSP = np.sort(np.abs(eigs))[[0, -1]]
        condKSP = lmaxKSP/lminKSP
        
        # Check convergence
        eu = errornorm(mms_data['u'], wh[0], 'H10')
        ez = errornorm(mms_data['flux'], wh[1], 'Hdiv')
        epT = errornorm(mms_data['pT'], wh[2], 'L2')
        ep = errornorm(mms_data['p'], wh[3], 'L2')

        h = mesh.hmin()
        ndofs = sum(Wi.dim() for Wi in W)
        errors = np.array([eu, ez, epT, ep])
        if errors0 is not None:
            rates = np.log(errors/errors0)/np.log(h/h0)
        else:
            rates = -np.ones_like(errors)
        errors0, h0 = errors, h

        ds_ = Measure('ds', domain=mesh, subdomain_data=boundaries)
        gd = sum(assemble(Constant(1)*ds_(tag)) for tag in u_dirichlet_tags)
        gp = sum(assemble(Constant(1)*ds_(tag)) for tag in p_dirichlet_tags)
        
        row = (h, ndofs, gd, gp) + sum(zip(errors, rates), ()) + (niters, lminKSP, lmaxKSP, condKSP)
        history.append(row)
        print(tabulate.tabulate(history, headers=headers))


        with open(get_path('cvrg', 'txt'), 'w') as out:
            np.savetxt(out, np.array(history), header=' '.join(headers))
