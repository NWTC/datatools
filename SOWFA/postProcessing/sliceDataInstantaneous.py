#
# Utilities to read in sliceDataInstantaneous legacy vtk output in ascii
# polydata format.
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import numpy as np

from vtk import vtkPolyDataReader
from vtk.util.numpy_support import vtk_to_numpy

def read_slice(fpath,const=None,verbose=True):
    """Read in vtk specified vtk file.

    If 'const' is not None, then check to make sure that the 'x','y',
    or 'z' cell-center values are constant. The returned arrays will
    have one less dimension.
    """
    dirmap = {'x':0, 'y':1, 'z':2}
    if const is not None:
        idir = dirmap[const.lower()]
    else:
        idir = -1
    reader = vtkPolyDataReader()
    reader.SetFileName(fpath)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()

    U = data.GetCellData().GetArray('U')
    U = vtk_to_numpy(U)
    Ncells = len(U)

    # get cell centers
    bounds = np.zeros((6,))
    cc = np.zeros(U.shape)
    for cell in range(Ncells):
        data.GetCellBounds(cell,bounds)
        cc[cell,0] = (bounds[0] + bounds[1])
        cc[cell,1] = (bounds[2] + bounds[3])
        cc[cell,2] = (bounds[4] + bounds[5])
    cc /= 2
    if verbose:
        for d,i in dirmap.items():
            print('{:s} in [{:g}, {:g}] m'.format(d,np.min(cc[:,i]),np.max(cc[:,i])))

    if const is not None:
        assert(np.max(cc[:,idir]) - np.min(cc[:,idir]) < 1e-8)

    # reorder (Fortran order, i.e., x increases fastest)
    order = np.arange(Ncells)
    sort0 = cc[:,0].argsort()
    cc = cc[sort0]
    sort1 = cc[:,1].argsort(kind='mergesort') # need stable sort so we don't lose the previous sorting
    cc = cc[sort1]
    sort2 = cc[:,2].argsort(kind='mergesort') # need stable sort so we don't lose the previous sorting
    cc = cc[sort2]
    order = order[sort0]
    order = order[sort1]
    order = order[sort2]
    U = U[order]

    # reshape arrays (assuming mesh is structured)
    xcc = cc[:,0]
    ycc = cc[:,1]
    zcc = cc[:,2]
    Nxy = len(np.nonzero(zcc==zcc[0])[0])
    Nxz = len(np.nonzero(ycc==ycc[0])[0])
    Nyz = len(np.nonzero(xcc==xcc[0])[0])
    assert(Nxy >= 1)
    assert(Nxz >= 1)
    assert(Nyz >= 1)
    Nx = np.sqrt(Nxz*Nxy/Nyz)
    Ny = Nxy / Nx
    Nz = Nyz / Ny
    assert(Nx==int(Nx))
    assert(Ny==int(Ny))
    assert(Nz==int(Nz))
    dims = [int(Nx),int(Ny),int(Nz)]
    if verbose:
        print('detected nx,ny,nz = {:d},{:d},{:d}'.format(*dims))

    if const is not None:
        assert(dims[idir] == 1)
        dims.pop(idir)
        if verbose:
            print('constant {:s} specified, new dimensions: {}'.format(const,dims))
    dims.append(3)
    x = cc.reshape(dims,order='F')
    U = U.reshape(dims,order='F')

    return x, U
