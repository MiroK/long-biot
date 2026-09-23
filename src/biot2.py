# NOTE: we assume the C param in (C + alpha^2/lmbda) to be C = 0 (which
# should be the most difficult case)

from collections import namedtuple
import argparse
from dolfin import *
from petsc4py import PETSc
import sympy as sp
import ulfy
from xii import *

from block.algebraic.petsc import KSP, LU, AMG
from block.block_mat import block_mat

print = PETSc.Sys.Print

from utils import StackOperator

from biot3 import (BiotParameters, Laplacian, generate_2d_domains, setup_2d_mms,
                   parse_V_bcs, parse_Q_bcs, get_path, SYM)

# ---

def get_system(boundaries, parameters, data, *, u_dirichlet_tags, p_dirichlet_tags, bdry_tags):
    '''Three field formulation'''
    mesh = boundaries.mesh()
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    n = FacetNormal(mesh)

    V = VectorFunctionSpace(mesh, 'CG', 2)
    Q = FunctionSpace(mesh, 'CG', 1)
    W = [V, Q]
    
    u, p = map(TrialFunction, W)
    v, q = map(TestFunction, W)

    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                              for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))

    a = block_form(W, 2)
    a[0][0] = inner(2*mu*SYM(grad(u)), SYM(grad(v)))*dx + inner(lmbda*div(u), div(v))*dx
    a[0][1] = -inner(alpha*p, div(v))*dx
    a[1][0] = -inner(alpha*q, div(u))*dx
    a[1][1] = -c*inner(p, q)*dx -inner(K*grad(p), grad(q))*dx

    u_neumann_tags = bdry_tags - set(u_dirichlet_tags)
    p_neumann_tags = bdry_tags - set(p_dirichlet_tags)

    L = block_form(W, 1)
    
    L[0] = inner(data['f_u'], v)*dx
    L[0] += sum(inner(data['u_neumann'][tag], v)*ds(tag) for tag in u_neumann_tags)
    
    # Pressure eq.
    L[1] = inner(-data['f_p'], q)*dx
    L[1] += sum(inner(data['p_neumann'][tag], q)*ds(tag) for tag in p_neumann_tags)

    # Nitsche bcs for pressure
    hF = CellDiameter(Q.mesh())
    nF = FacetNormal(Q.mesh())
    gammaF = Constant(5)
    for tag in p_dirichlet_tags:
        a[1][1] += (+ inner(dot(K*grad(p), nF), q)*ds(tag)
                    + inner(dot(K*grad(q), nF), p)*ds(tag)
                    - (K*gammaF/hF)*inner(p, q)*ds(tag))

        p0 = data['p_dirichlet'][tag]
        L[1] += (+ inner(dot(K*grad(q), nF), p0)*ds(tag)
                 - (K*gammaF/hF)*inner(p0, q)*ds(tag))

    V_bcs = [DirichletBC(V, data['u_dirichlet'][tag], boundaries, tag) for tag in u_dirichlet_tags]
    Q_bcs = []
    W_bcs = [V_bcs, Q_bcs]

    A, b = map(ii_assemble, (a, L))
    A, b = apply_bc(A, b, bcs=W_bcs)

    return A, b, W, W_bcs


def get_inner_product_standard(boundaries, parameters, *, u_dirichlet_tags, p_dirichlet_tags, W, Wbcs, bdry_tags,
                               inverseQ='lu'):
    '''Structure V x (Q x QT) where pressures are coupled'''
    V, Q = W
    u, p = map(TrialFunction, W)
    v, q = map(TestFunction, W)

    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                               for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))


    b0 = inner(2*mu*SYM(grad(u)), SYM(grad(v)))*dx + (1+lmbda)*inner(div(u), div(v))*dx
    L = inner(Constant((0, 0)), v)*dx
    B0, _ = assemble_system(b0, L, Wbcs[0])

    ds = Measure('ds', domain=Q.mesh(), subdomain_data=boundaries)
    # ----
    b1 = (c + alpha**2/(1+lmbda))*inner(p, q)*dx + inner(K*grad(p), grad(q))*dx
    # Nietsche terms
    hF = CellDiameter(Q.mesh())
    nF = FacetNormal(Q.mesh())
    gammaF = Constant(5)
    for tag in p_dirichlet_tags:
        b1 += (- inner(dot(K*grad(p), nF), q)*ds(tag)
               - inner(dot(K*grad(q), nF), p)*ds(tag)
               + (K*gammaF/hF)*inner(p, q)*ds(tag))    
    

    B1 = ii_assemble(b1)

    BB = block_diag_mat([B0, B1])
    # In inversion we do 1 x 2
    iB00 = LU(B0)

    if inverseQ == 'lu':
        iB11 = LU(B1)
    elif inverseQ == 'amg':
        raise ValueError
        iB11 = AMG(B1)
    iBB = block_diag_mat([iB00, iB11])

    return BB, iBB


def get_inner_product_espen(boundaries, parameters, *, u_dirichlet_tags, p_dirichlet_tags, W, Wbcs, bdry_tags,
                                     inverseQ='lu'):
    '''Structure V x (Q x QT) where pressures are coupled'''
    V, Q = W
    u, p = map(TrialFunction, W)
    v, q = map(TestFunction, W)

    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                               for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))

    
    b0 = inner(2*mu*SYM(grad(u)), SYM(grad(v)))*dx + (1+lmbda)*inner(div(u), div(v))*dx
    L = inner(Constant((0, 0)), v)*dx
    B0, _ = assemble_system(b0, L, Wbcs[0])


    ds = Measure('ds', domain=Q.mesh(), subdomain_data=boundaries)        
    hF = CellDiameter(Q.mesh())
    nF = FacetNormal(Q.mesh())
    gammaF = Constant(5)
    
    def c_form(p, q, ds=ds, hF=hF, gammaF=gammaF):
        # Now components of Espen preconditioner
        a = c*inner(p, q)*dx + inner(K*grad(p), grad(q))*dx

        # Pressure bcs
        for tag in p_dirichlet_tags:
            a += (- inner(dot(K*grad(p), nF), q)*ds(tag)
                  - inner(dot(K*grad(q), nF), p)*ds(tag)
                  + (K*gammaF/hF)*inner(p, q)*ds(tag))    
        return a

    # Extended space for the inverse
    QQ = FunctionSpace(mesh, MixedElement([Q.ufl_element()]*2))

    p0, p1 = TrialFunctions(QQ)
    q0, q1 = TestFunctions(QQ)
    
    # ---

    bc_tags = {'T': set(bdry_tags) - set(u_dirichlet_tags),
               'P': set()}
    scale = Constant(1) # FIXME, this will be the thickness
    kappa = alpha**2/(1+lmbda)*scale**2
    k_form, ker = Laplacian((p1, q1), boundaries, bc_tags, kappa=kappa)

    # ---
    m_form = (alpha**2/(1+lmbda))*inner(p0, q0)*dx

    e_form = c_form(p0, q0) + c_form(p1, q0) + c_form(p0, q1) + c_form(p1, q1)
    e_form += m_form + k_form

    EE = assemble(e_form)

    # ----

    precond0 = LU(B0)

    if inverseQ == 'lu':
        invE = LU(EE)
    elif inverseQ == 'amg':

        invE = AMG(EE,
                   parameters={
                       'pc_hypre_boomeramg_strong_threshold': 0.1,
                       'pc_hypre_boomeramg_nodal_coarsen': 1,
                       'pc_hypre_boomeramg_vec_interp_variant': 1,
                       'pc_hypre_boomeramg_interp_type': 'ext+i',
                       'pc_hypre_boomeramg_smooth_type': 'Schwarz-smoothers'
                   })

    S = StackOperator(Q, QQ)  

    precond1 = S.T*invE*S
    
    iBB = block_diag_mat([precond0, precond1])

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
    parser.add_argument('-inverseQ', type=str, default='lu', choices=('lu', 'amg'))
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


    result_dir = f'./results/biot2/precond{args.precond}'
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
               '|eu|1', 'reu', '|ep|1', 'rep1',
               'niters', 'lminKSP', 'lmaxKSP', 'condKSP')

    length = args.L
    
    h0, errors0, history = None, None, []
    for (level, boundaries) in zip(range(args.nrefs), mesh_gen):
        mesh = boundaries.mesh()
        
        A, b, W, W_bcs = get_system(boundaries, parameters=parameters, data=mms_data,
                                    u_dirichlet_tags=u_dirichlet_tags,
                                    p_dirichlet_tags=p_dirichlet_tags,
                                    bdry_tags=boundary_tags)

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
        ep = errornorm(mms_data['p'], wh[1], 'H10')

        h = mesh.hmin()
        ndofs = sum(Wi.dim() for Wi in W)
        errors = np.array([eu, ep])
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
