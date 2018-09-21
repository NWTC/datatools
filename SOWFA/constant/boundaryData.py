#
# Module to handle SOWFA boundary data that belong in
#   casedir/constant/boundaryData
#
# Written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import os
import numpy as np


pointsheader = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.x                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       vectorField;
    location    "constant/boundaryData/{patchName}";
    object      points;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{N}
("""

dataheader = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.x                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {patchType}AverageField;
    location    "constant/boundaryData/{patchName}/{timeName}";
    object      values;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Average
{avgValue}

{N}
("""


def write_points(fname,x,y,z,patchName='patch'):
    """Write out a points file which should be stored in
        constant/boundaryData/patchName/points
    """
    N = len(x)
    assert(N == len(y) == len(z))
    if len(x.shape) > 1:
        x = x.ravel()
        y = y.ravel()
        z = z.ravel()
        N = len(x)
    dpath = os.path.split(fname)[0]
    if not os.path.isdir(dpath):
        os.makedirs(dpath)
    np.savetxt(fname,
               np.stack((x,y,z)).T, fmt='(%f %f %f)',
               header=pointsheader.format(patchName=patchName,N=N),
               footer=')',
               comments='')

def write_data(fname,
               data,
               patchName='patch',
               timeName=0,
               avgValue=None):
    """Write out a boundaryData file which should be stored in
        constant/boundarydata/patchName/timeName/fieldName

    Parameters
    ----------
    fname : str
        Output data file name
    data : numpy.ndarray
        Field data to be written out, with shape (3,N) for vectors
        and shape (N) for scalars; 2-D or 3-D data should be flattened
        beforehand.
    patchName : str, optional
        Name of the boundary patch
    timeName : scalar or str, optional
        Name of corresponding time directory
    avgValue : scalar or list-like, optional
        To set avgValue in the data file; probably not used.

    @author: ewquon
    """
    dims = data.shape
    N = dims[-1]
    if len(dims) == 1:
        patchType = 'scalar'
        if avgValue is None:
            avgValueStr = '0'
        else:
            avgValueStr = str(avgValue)
    elif len(dims) == 2:
        patchType = 'vector'
        assert(dims[0] == 3)
        if avgValue is None:
            avgValueStr = '(0 0 0)'
        else:
            avgValueStr = '(' + ' '.join([str(val) for val in list(avgValue)]) + ')'
    else:
        print('ERROR: Unexpected number of dimensions! No data written.')
        return

    dpath = os.path.split(fname)[0]
    if not os.path.isdir(dpath):
        os.makedirs(dpath)

    headerstr = dataheader.format(patchType=patchType,
                                  patchName=patchName,
                                  timeName=timeName,
                                  avgValue=avgValueStr,
                                  N=N)
    if patchType == 'vector':
        np.savetxt(fname,
                   data.T, fmt='(%g %g %g)',
                   header=headerstr, footer=')',
                   comments='')
    elif patchType == 'scalar':
        np.savetxt(fname,
                   data.reshape((N,1)), fmt='%g',
                   header=headerstr, footer=')',
                   comments='')


class CartesianPatch(object):
    """Object to facilitate outputing boundary patches on a Cartesian
    grid, with boundaries at x=const or y=const.
    """

    def __init__(self,x,y,z,dpath='.',name='patch'):
        """For a Cartesian mesh, the grid is defined by 1-D coordinate
        vectors x,y,z
        """
        self.dpath = dpath
        self.name = name
        # check for constant x/y
        if np.all(x==x[0]):
            self.desc = 'const x={:.1f}'.format(x[0])
            x = [x[0]]
        elif np.all(y==y[0]):
            self.desc = 'const y={:.1f}'.format(y[0])
            y = [y[0]]
        else:
            raise ValueError('x and y not constant, domain is not Cartesian')
        # set up points
        self.x = x
        self.y = y
        self.z = z
        self.Nx = len(x)
        self.Ny = len(y)
        self.Nz = len(z)
        # set up mesh
        self.X, self.Y, self.Z = np.meshgrid(x,y,z,indexing='ij')

    def __repr__(self):
        s = 'Cartesian patch "{:s}" : '.format(self.name)
        s += self.desc
        s +=', size=({:d},{:d},{:d})'.format(self.Nx,self.Ny,self.Nz)
        return s

    def write_points(self,fpath=None):
        """Write out constant/boundaryData/patchName/points file"""
        if fpath is None:
            fpath = os.path.join(self.dpath,self.name,'points')
        write_points(fpath,self.X,self.Y,self.Z,patchName=self.name)
        print('Wrote points to '+fpath)

    def write_profiles(self, t, z,
                       U=None, V=None, W=None, T=None, k=None,
                       time_range=(None,None),
                       verbose=True):
        """Write out constant/boundaryData/patchName/*/{U,T} given a
        set of time-height profiles. Outputs will be interpolated to 
        the patch heights.

        Inputs
        ------
        t : np.ndarray
            Time vector with length Nt
        z : np.ndarray
            Height vector with length Nz
        U : np.ndarray
            Velocity vectors (Nt,Nz,3) or x-velocity component (Nt,Nz)
        V : np.ndarray, optional
            y-velocity component (Nt,Nz)
        W : np.ndarray, optional
            z-velocity component (Nt,Nz)
        T : np.ndarray
            potential temperature (Nt,Nz)
        time_range : tuple, optional
            range of times to write out boundary data
        """
        assert(U is not None)
        assert(T is not None)
        Nt, Nz = U.shape[:2]
        assert(Nt == len(t))
        assert(Nz == len(z))
        if len(U.shape) == 2:
            have_components = True
            assert(np.all(U.shape == (Nt,Nz)) and \
                    np.all(V.shape == (Nt,Nz)) and \
                    np.all(W.shape == (Nt,Nz)) and \
                    np.all(T.shape == (Nt,Nz)))
        elif len(U.shape) == 3:
            have_components = False
            assert(U.shape[2] == 3)
            assert(np.all(T.shape == (Nt,Nz)))
        else:
            raise InputError

        if not have_components:
            V = U[:,:,1]
            W = U[:,:,2]
            U = U[:,:,0]

        if k is not None:
            assert(np.all(k.shape == (Nt,Nz)))

        if time_range[0] is None:
            time_range[0] = 0.0
        if time_range[1] is None:
            time_range[1] = 9e9

        if (not len(self.z) == Nz) or (not np.all(self.z == z)):
            interpolate = True
        else:
            # input z are equal to patch z
            interpolate = False

        Upatch = np.zeros((self.Nx,self.Ny,self.Nz))
        Vpatch = np.zeros((self.Nx,self.Ny,self.Nz))
        Wpatch = np.zeros((self.Nx,self.Ny,self.Nz))
        Tpatch = np.zeros((self.Nx,self.Ny,self.Nz))
        kpatch = np.zeros((self.Nx,self.Ny,self.Nz))
        for it,ti in enumerate(t):
            if (ti < time_range[0]) or (ti > time_range[1]):
                if verbose:
                    print('skipping time '+str(ti))
                continue
            tname = '{:f}'.format(ti).rstrip('0').rstrip('.')

            timepath = os.path.join(self.dpath, self.name, tname)
            if not os.path.isdir(timepath):
                os.makedirs(timepath)

            Upatch[:,:,:] = 0.0
            Vpatch[:,:,:] = 0.0
            Wpatch[:,:,:] = 0.0
            Tpatch[:,:,:] = 0.0
            kpatch[:,:,:] = 0.0
            if interpolate:
                if verbose:
                    print('interpolating data at t = {:s} s'.format(tname)) 
                for iz,zp in enumerate(self.z):
                    i = np.nonzero(z > zp)[0][0]
                    f = (zp - z[i-1]) / (z[i] - z[i-1])
                    Upatch[:,:,iz] = U[it,i-1] + f*(U[it,i]-U[it,i-1])
                    Vpatch[:,:,iz] = V[it,i-1] + f*(V[it,i]-V[it,i-1])
                    Wpatch[:,:,iz] = W[it,i-1] + f*(W[it,i]-W[it,i-1])
                    Tpatch[:,:,iz] = T[it,i-1] + f*(T[it,i]-T[it,i-1])
                    if k is not None:
                        kpatch[:,:,iz] = k[it,i-1] + f*(k[it,i]-k[it,i-1])
            else:
                if verbose:
                    print('mapping data at t = {:s} s'.format(tname)) 
                for iz in range(Nz):
                    Upatch[:,:,iz] = U[it,iz]
                    Vpatch[:,:,iz] = V[it,iz]
                    Wpatch[:,:,iz] = W[it,iz]
                    Tpatch[:,:,iz] = T[it,iz]
                    if k is not None:
                        kpatch[:,:,iz] = k[it,iz]
        
            Upath = os.path.join(timepath,'U')
            Udata = np.stack((Upatch.ravel(), Vpatch.ravel(), Wpatch.ravel()))
            if verbose: print('writing data to {:s}'.format(Upath)) 
            write_data(Upath, Udata, patchName=self.name, timeName=ti)

            Tpath = os.path.join(timepath,'T')
            if verbose: print('writing data to {:s}'.format(Tpath)) 
            write_data(Tpath, Tpatch.ravel(), patchName=self.name, timeName=ti)

            if k is not None:
                kpath = os.path.join(timepath,'k')
                if verbose: print('writing data to {:s}'.format(kpath)) 
                write_data(kpath, kpatch.ravel(), patchName=self.name, timeName=ti)


