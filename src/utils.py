from block.block_base import block_base
from block.object_pool import vec_pool
from block import block_vec

from dolfin import FunctionSpace, MixedElement, Function

from scipy.sparse import csr_matrix
from xii.linalg.convert import numpy_to_petsc
import numpy as np


def get_domain_diameter(mesh):
    '''Well if it is rectangle'''
    xmin, ymin = mesh.coordinates().min(axis=0)
    xmax, ymax = mesh.coordinates().max(axis=0)
    return max(ymax - ymin, xmax - xmin)


def StackOperator(dspaces, rspaces):
    '''Q -> Q x Q'''
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


class SerializeOperator(block_base):
    '''Reordering'''
    def __init__(self, mspace):
        assert isinstance(mspace, FunctionSpace)

        elm = mspace.ufl_element()
        assert isinstance(elm, MixedElement)
        assert elm.num_sub_elements() > 0

        self.dofs = [mspace.sub(i).dofmap().dofs() for i in range(elm.num_sub_elements())]
        self.W = mspace

    @vec_pool
    def create_vec(self, dim=1):
        return Function(self.W).vector()

    def matvec(self, b):
        x = self.create_vec(dim=0)
        x_arr = x.get_local()

        start = 0
        b_arr = b.get_local()
        for (xi, dofsi) in zip(x, self.dofs):
            x_arr[start:start+len(dofsi)] = b_arr[dofsi]
            start += len(dofsi)
        x.set_local(x_arr)
        return x

    def transpmult(self, b):
        x = self.create_vec(dim=0)
        x_arr = x.get_local()

        start = 0
        b_arr = b.get_local()
        for (xi, dofsi) in zip(x, self.dofs):
            x_arr[dofsi] = b_arr[start:start+len(dofsi)]
            start += len(dofsi)
        x.set_local(x_arr)
        return x

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
