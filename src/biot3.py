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
import os

print = PETSc.Sys.Print

BiotParameters = namedtuple('BiotParameters', ('alpha', 'K', 'mu', 'lmbda', 'c'))


SYM = lambda x: sym(x)


def Laplacian(arg, boundaries, bc_tags, kappa):
    '''-div(kappa*grad) with Dirichlet boundaries on tagged parts'''
    if isinstance(arg, FunctionSpace):
        V = arg
        u, v = TrialFunction(V), TestFunction(V)
        mesh = V.mesh()        
    else:
        u, v = arg
        mesh = u.ufl_domain().ufl_cargo()

    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    dx = Measure('dx', domain=mesh)
    
    hFi, hFe = CellDiameter(mesh), CellDiameter(mesh)

    a = kappa*inner(grad(u), grad(v))*dx + kappa*(1/avg(hFi))*inner(jump(u), jump(v))*dS

    volume = assemble(Constant(1)*dx)
    
    ker = [Constant(1/sqrt(volume))]
    if bc_tags['P'] or bc_tags['T']:
        ker.pop()

    # We have bcs and no kernel
    for tag in (bc_tags['P'] | bc_tags['T']):
        print(tag)
        a += kappa*(1/hFi)*inner(u, v)*ds(tag)

    return (a, ker)



def generate_2d_domains(L, ncells0, kind='crossed'):
    '''Unit square with marked edges'''
    subdomains = {1: CompiledSubDomain('near(x[0], 0)'),
                  2: CompiledSubDomain('near(x[0], L)', L=L),
                  3: CompiledSubDomain('near(x[1], 0)'),
                  4: CompiledSubDomain('near(x[1], 1)')}

    boundary_tags = set((1, 2, 3, 4))

    yield boundary_tags
    
    ncells = ncells0
    while True:
        mesh = RectangleMesh(Point(0, 0), Point(L, 1), L*ncells, ncells, kind)

        boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
        for (tag, subd) in subdomains.items():
            subd.mark(boundaries, tag)

        yield boundaries
            
        ncells *= 2
    

def setup_2d_mms(parameters):
    '''Manufactured problem on (0, 1)^2'''
    x, y = SpatialCoordinate(UnitSquareMesh(2, 2))

    u = as_vector((sin(pi*(x+y)), sin(2*pi*(x-y))))
    p = cos(pi*(x-y))

    mu, K, alpha, lmbda = Constant(1), Constant(1), Constant(1), Constant(1)
    c = Constant(1)
    
    pT = -lmbda*div(u) + alpha*p
    stress = 2*mu*SYM(grad(u)) - pT*Identity(2)
    flux = -K*grad(p)

    f_u = -div(stress)
    f_p = c*p + alpha*div(u) + div(flux)

    normals = [Constant((-1, 0)), Constant((1, 0)), Constant((0, -1)), Constant((0, 1))]

    stress_data = [dot(stress, n) for n in normals]
    flux_data = [dot(flux, n) for n in normals]

    mu0, K0, alpha0, lmbda0, c0 = sp.symbols('mu K alpha lmbda c')
    subs = {mu: mu0, K: K0, alpha: alpha0, lmbda: lmbda0, c: c0}

    def as_expr(v, subs=subs, params=parameters):
        expr = ulfy.Expression(v, degree=5, subs=subs)
        for key in expr._user_parameters:
            if hasattr(params, key):
                setattr(expr, key, getattr(params, key))
        return expr

    return {'u': as_expr(u), 'p': as_expr(p), 'pT': as_expr(pT),
            'f_u': as_expr(f_u),
            'f_p': as_expr(f_p),
            'u_neumann': dict(enumerate(map(as_expr, stress_data), 1)),
            'p_neumann': dict(enumerate(map(as_expr, flux_data), 1)),
            'u_dirichlet': dict(enumerate(map(as_expr, [u]*4), 1)),
            'p_dirichlet': dict(enumerate(map(as_expr, [p]*4), 1))}


# ---

def get_system(boundaries, parameters, data, *, u_dirichlet_tags, p_dirichlet_tags, bdry_tags,
               pdegrees='2_2_1'):
    '''Three field formulation'''
    mesh = boundaries.mesh()
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    n = FacetNormal(mesh)

    V_deg, Q_deg, QT_deg = map(int, pdegrees.split('_'))
    V = VectorFunctionSpace(mesh, 'CG', V_deg)
    Q = FunctionSpace(mesh, 'CG', Q_deg)
    QT = FunctionSpace(mesh, 'CG', QT_deg)
    W = [V, Q, QT]
    
    u, p, pT = map(TrialFunction, W)
    v, q, qT = map(TestFunction, W)

    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                              for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))

    a = block_form(W, 2)
    a[0][0] = inner(2*mu*SYM(grad(u)), SYM(grad(v)))*dx
    a[0][2] = -inner(pT, div(v))*dx
    a[1][1] = -(c + alpha**2/lmbda)*inner(p, q)*dx -inner(K*grad(p), grad(q))*dx
    a[1][2] = + (alpha/lmbda)*inner(pT, q)*dx
    a[2][0] = -inner(qT, div(u))*dx
    a[2][1] = + (alpha/lmbda)*inner(p, qT)*dx
    a[2][2] = -(1/lmbda)*inner(pT, qT)*dx

    u_neumann_tags = bdry_tags - set(u_dirichlet_tags)
    p_neumann_tags = bdry_tags - set(p_dirichlet_tags)

    L = block_form(W, 1)
    
    L[0] = inner(data['f_u'], v)*dx
    L[0] += sum(inner(data['u_neumann'][tag], v)*ds(tag) for tag in u_neumann_tags)
    
    # Pressure eq.
    L[1] = inner(-data['f_p'], q)*dx
    L[1] += sum(inner(data['p_neumann'][tag], q)*ds(tag) for tag in p_neumann_tags)

    V_bcs = [DirichletBC(V, data['u_dirichlet'][tag], boundaries, tag) for tag in u_dirichlet_tags]

    # Nietsche bcs
    nF, hF = FacetNormal(mesh), CellDiameter(mesh)
    gammaF = Constant(5)
    for tag in p_dirichlet_tags:
        a[1][1] += (
            inner(dot(K*grad(p), nF), q)*ds(tag)
            + inner(dot(K*grad(q), nF), p)*ds(tag)
            - (K*gammaF/hF)*inner(p, q)*ds(tag)
        )

        p0 = data['p_dirichlet'][tag]

        L[1] += (
            inner(dot(K*grad(q), nF), p0)*ds(tag)
            - (K*gammaF/hF)*inner(p0, q)*ds(tag)
        )
    
    Q_bcs = []
    QT_bcs = []
    W_bcs = [V_bcs, Q_bcs, QT_bcs]

    A, b = map(ii_assemble, (a, L))
    A, b = apply_bc(A, b, bcs=W_bcs)
    
    return A, b, W, W_bcs


def get_inner_product_standard(boundaries, parameters, *, u_dirichlet_tags, p_dirichlet_tags, W, Wbcs, bdry_tags):
    '''Structure V x (Q x QT) where pressures are coupled'''
    u, p, pT = map(TrialFunction, W)
    v, q, qT = map(TestFunction, W)

    V, Q, QT = W

    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                               for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))

    a = block_form(W, 2)
    a[0][0] = inner(2*mu*SYM(grad(u)), SYM(grad(v)))*dx
    a[1][1] = (c + alpha**2/lmbda)*inner(p, q)*dx + inner(K*grad(p), grad(q))*dx

    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    nF, hF = FacetNormal(Q.mesh()), CellDiameter(Q.mesh())
    gammaF = Constant(5)
    # Add Nietsche terms
    for tag in p_dirichlet_tags:
        a[1][1] += (
            - inner(dot(K*grad(p), nF), q)*ds(tag)
            - inner(dot(K*grad(q), nF), p)*ds(tag)
            + (K*gammaF/hF)*inner(p, q)*ds(tag)
        )
    
    a[1][2] = - (alpha/lmbda)*inner(pT, q)*dx
    a[2][1] = - (alpha/lmbda)*inner(p, qT)*dx
    a[2][2] = (1/2/mu + 1/lmbda)*inner(pT, qT)*dx

    B = ii_assemble(a)
    B, _ = apply_bc(B, b=None, bcs=[Wbcs[0], [], []])

    # In inversion we do 1 x 2
    B00 = B[0][0]
    B11 = block_mat([[B[1][1], B[1][2]],
                     [B[2][1], B[2][2]]])

    iB00 = LU(B00)
    iB11 = LU(monolithic(B11))
    R = ReductionOperator([1, 3], W)
    iBB = R.T*block_diag_mat([iB00, iB11])*R

    return B, iBB


def get_inner_product_espen(boundaries, parameters, *, u_dirichlet_tags, p_dirichlet_tags, W, Wbcs, bdry_tags):
    '''Structure V x (Q x QT) where pressures are coupled'''
    u, p, pT = map(TrialFunction, W)
    v, q, qT = map(TestFunction, W)

    V, Q, QT = W
    
    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                               for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))

    bV = inner(2*mu*SYM(grad(u)), SYM(grad(v)))*dx
    lV = inner(Constant((0, 0)), v)*dx
    B0, _ = assemble_system(bV, lV, Wbcs[0])

    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

    QQ = [Q, QT, QT]
    # ---
    a = block_form(QQ, 2)
    a[0][0] = (c + alpha**2/lmbda)*inner(p, q)*dx + inner(K*grad(p), grad(q))*dx

    nF, hF = FacetNormal(Q.mesh()), CellDiameter(Q.mesh())
    gammaF = Constant(5)
    # Add Nietsche terms
    for tag in p_dirichlet_tags:
        a[0][0] += (
            - inner(dot(K*grad(p), nF), q)*ds(tag)
            - inner(dot(K*grad(q), nF), p)*ds(tag)
            + (K*gammaF/hF)*inner(p, q)*ds(tag)
        )
        
    
    a[0][1] = (-alpha/lmbda)*inner(pT, q)*dx
    a[0][2] = (-alpha/lmbda)*inner(pT, q)*dx

    a[1][0] = (-alpha/lmbda)*inner(qT, p)*dx

    a[1][1] = (1/2/mu + 1/lmbda)*inner(pT, qT)*dx
    a[1][2] = (1/lmbda)*inner(pT, qT)*dx

    a[2][0] = (-alpha/lmbda)*inner(qT, p)*dx    
    a[2][1] = (1/lmbda)*inner(pT, qT)*dx

    bc_tags = {'T': set(bdry_tags) - set(u_dirichlet_tags),
               'P': set()}
    scale = Constant(1) # FIXME, this will be the thickness
    kappa = scale**2/2/mu
    k_form, ker = Laplacian(QT, boundaries, bc_tags, kappa=kappa)

    a[2][2] = k_form + (1/lmbda)*inner(pT, qT)*dx

    precond0 = LU(B0)
    
    # Puttin together
    EE = ii_assemble(a)

    E = monolithic(EE)
    invE = LU(E)
    precond1 = invE

    S = StackOperator(2, QT, W=[V, Q])
    R = ReductionOperator([1, 4], [V, Q, QT, QT])

    precond = block_diag_mat([precond0, precond1])
    
    iBB = S.T*R.T*precond*R*S

    return None, iBB


def get_inner_product_espen_diagonal(boundaries, parameters, *, u_dirichlet_tags, p_dirichlet_tags, W, Wbcs, bdry_tags):
    '''Structure V x (Q x QT) where pressures are coupled'''
    u, p, pT = map(TrialFunction, W)
    v, q, qT = map(TestFunction, W)

    V, Q, QT = W
    
    mu, lmbda, alpha, K, c = (Constant(getattr(parameters, p))
                               for p in ('mu', 'lmbda', 'alpha', 'K', 'c'))

    bV = inner(2*mu*SYM(grad(u)), SYM(grad(v)))*dx
    lV = inner(Constant((0, 0)), v)*dx
    B0, _ = assemble_system(bV, lV, Wbcs[0])

    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

    # ---

    b1 = (c + alpha**2/lmbda)*inner(p, q)*dx + inner(K*grad(p), grad(q))*dx

    nF, hF = FacetNormal(Q.mesh()), CellDiameter(Q.mesh())
    gammaF = Constant(5)
    # Add Nietsche terms
    for tag in p_dirichlet_tags:
        a[0][0] += (
            - inner(dot(K*grad(p), nF), q)*ds(tag)
            - inner(dot(K*grad(q), nF), p)*ds(tag)
            + (K*gammaF/hF)*inner(p, q)*ds(tag)
        )
    B1 = assemble(b1)
    

    QQ = [QT, QT]
    
    a = block_form(QQ, 2)    
    a[0][0] = (1/2/mu + 1/lmbda)*inner(pT, qT)*dx
    a[0][1] = (1/lmbda)*inner(pT, qT)*dx

    a[1][0] = (1/lmbda)*inner(pT, qT)*dx

    bc_tags = {'T': set(bdry_tags) - set(u_dirichlet_tags),
               'P': set()}
    scale = Constant(1) # FIXME, this will be the thickness
    kappa = scale**2/2/mu
    k_form, ker = Laplacian(QT, boundaries, bc_tags, kappa=kappa)

    a[1][1] = k_form + (1/lmbda)*inner(pT, qT)*dx

    EE = ii_assemble(a)
    
    # Puttin together
    precond0 = LU(B0)
    precond1 = LU(B1)

    E = monolithic(EE)
    invE = LU(E)
    precond2 = invE

    S = StackOperator(2, QT, W=[V, Q])
    R = ReductionOperator([1, 2, 4], [V, Q, QT, QT])

    precond = block_diag_mat([precond0, precond1, precond2])
    
    iBB = S.T*R.T*precond*R*S

    return None, iBB

# ---

def parse_V_bcs(bcs, tags):
    '''Only return displacement'''
    assert len(bcs) == len(tags)
    assert all(bc in ('D', 'T') for bc in bcs)

    return [tag for (tag, bc) in zip(tags, bcs) if bc == 'D']


def parse_Q_bcs(bcs, tags):
    '''Only return displacement'''
    assert len(bcs) == len(tags)
    assert all(bc in ('P', 'F') for bc in bcs)

    return [tag for (tag, bc) in zip(tags, bcs) if bc == 'P']


def get_path(what, ext, result_dir, ignore_keys, args):
    template_path = '_'.join([what] + [f'{key.upper()}{value}' for key, value in args.items()
                                       if key not in ignore_keys])
    template_path = '.'.join([template_path, ext])
    path = os.path.join(result_dir, template_path)
    print(f'Saving to {path}')
    return path

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


    result_dir = f'./results/biot3/precond{args.precond}'
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


    get_inner_product = {'standard': get_inner_product_standard,
                         'espen': get_inner_product_espen,
                         'espendiag': get_inner_product_espen_diagonal}[args.precond]

    parameters = BiotParameters(alpha=args.alpha, K=args.K, mu=args.mu, lmbda=args.lmbda, c=args.c)
    
    mms_data = setup_2d_mms(parameters)

    # Common solver settings
    opts = PETSc.Options()        
    opts.setValue('ksp_rtol', 1E-12)
    opts.setValue('ksp_view_pre', None)
    opts.setValue('ksp_monitor_true_residual', None)
    opts.setValue('ksp_converged_reason', None)

    headers = ('h', 'ndofs', '|GD|', '|GP|',
               '|eu|1', 'reu', '|ep|1', 'rep1', '|epT|0', 'repT',
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
                                    bdry_tags=boundary_tags)

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
        epT = errornorm(mms_data['pT'], wh[2], 'L2')

        h = mesh.hmin()
        ndofs = sum(Wi.dim() for Wi in W)
        errors = np.array([eu, ep, epT])
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
