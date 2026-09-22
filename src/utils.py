from block import block_mat
from dolfin import FunctionSpace

from scipy.sparse import csr_matrix
from xii.linalg.convert import numpy_to_petsc
import numpy as np


def StackOperator(dspaces, rspaces):
    if isinstance(dspaces, FunctionSpace):
        assert isinstance(rspaces, FunctionSpace)

        delm = dspaces.ufl_element()

        relm, = set(rspaces.ufl_element().sub_elements())
        assert delm == relm

        nblocks = rspaces.ufl_element().num_sub_elements()
        assert nblocks > 0
        
        rows = np.arange(rspaces.dim())
        cols = np.repeat(np.arange(dspaces.dim()), nblocks)
        data = np.ones_like(cols)

        S = csr_matrix((data, (rows, cols)), shape=(rspaces.dim(), dspaces.dim()))

        return numpy_to_petsc(S)

    raise ValueError
        
# -------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *


    mesh = UnitSquareMesh(4, 4)

    elm = FiniteElement('Lagrange', mesh.ufl_cell(), 1)

    V = FunctionSpace(mesh, elm)
    VV = FunctionSpace(mesh, MixedElement([elm]*3))


    S = StackOperator(V, VV)

    fV = interpolate(Expression('x[0]+x[1]', degree=1), V)

    fVV = Function(VV)
    fVV.vector()[:] = S*fV.vector()

    fVis = split(fVV)

    for fVi in fVis:
        c = inner(fV, fV)*dx
        e = inner(fV - fVi, fV - fVi)*dx
        print(sqrt(abs(assemble(e))), sqrt(abs(assemble(c))))
    
