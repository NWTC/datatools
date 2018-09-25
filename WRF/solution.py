#
# Process WRF solution file
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import sys, os
import numpy as np

#from netCDF4 import Dataset
from netCDF4 import Dataset, MFDataset
try:
   import xarray
except ImportError:
    have_xarray = False
else:
    print('xarray reader available')
    have_xarray = True

g = 9.81
default_aggdim = 'time'

class WRFSolution(object):
    """Object to hold a single WRF solution snapshot"""

    def __init__(self,*args,**kwargs):
        verbose = kwargs.get('verbose',True)
        aggdim = kwargs.get('aggdim',default_aggdim)
        self.use_xarray = kwargs.get('use_xarray',have_xarray)
        if self.use_xarray:
            desc = 'with xarray'
        else:
            desc = 'with netcdf'
        Nfiles = len(args)
        self.filelist = []
        for fpath in [ fpath for fpath in args if os.path.isfile(fpath) ]:
            try:
                Dataset(fpath)
            except (IOError,OSError): # NetCDF: Unknown file format
                pass
            else:
                self.filelist.append(fpath)
        if self.use_xarray:
            nc = xarray.open_mfdataset(self.filelist, concat_dim=aggdim)
            self.Nt, self.Nz, self.Ny, self.Nx = nc.variables['U'].shape
            self.Nx -= 1 # U is staggered in x
        else:
            nc = MFDataset(self.filelist, aggdim=aggdim)
            self.Nt = len(nc.dimensions['time'])
            self.Nx = len(nc.dimensions['west_east'])
            self.Ny = len(nc.dimensions['south_north'])
            self.Nz = len(nc.dimensions['bottom_top'])
        self.varlist = list(nc.variables)
        self._read_vars(nc)

    def _read_vars(self,nc):
        # unstaggered
        self.T = nc.variables['T'][:] + 300.0
        # staggered in x
        U = nc.variables['U'][:]
        self.U = 0.5*(U[:,:,:,:-1] + U[:,:,:,1:])
        # staggered in y
        V = nc.variables['V'][:]
        self.V = 0.5*(V[:,:,:-1,:] + V[:,:,1:,:])
        # staggered in z
        W = nc.variables['W'][:]
        PH = nc.variables['PH'][:]
        PHB = nc.variables['PHB'][:]
        self.W = 0.5*(W[:,:-1,:,:] + W[:,1:,:,:])
        # calculate z == (ph + phb)/g
        self.z = 0.5*( PH[:,:-1,:,:] +  PH[:,1:,:,:] +
                      PHB[:,:-1,:,:] + PHB[:,1:,:,:] ) / g

        # calculate height AGL
        if 'HGT' in self.varlist:
            # TODO: test this
            hgt = nc.variables['HGT'][:]
            for i in range(self.Nx):
                for j in range(self.Ny):
                    self.z[:,i,j,:] -= hgt[i,j]

        # xarray doesn't read in the mfdataset until we call .values
        if self.use_xarray:
            self.z = self.z.values
            self.U = self.U.values
            self.V = self.V.values
            self.W = self.W.values
            self.T = self.T.values

    def sample_profile(self,itime=slice(0,None),i=None,j=None,overwrite=False):
        """Extracts velocity and temperature profile at a specified
        location (defaults to center of domain).
        
        If overwrite is True, reduce the dimensions of the stored z, U,
        V, and T variables; otherwise, return the profiles.
        """
        if i is None:
            i = int(self.Nx / 2)
        if j is None:
            j = int(self.Ny / 2)
        zprofile = self.z[itime,:,j,i]
        Uprofile = self.U[itime,:,j,i]
        Vprofile = self.V[itime,:,j,i]
        Wprofile = self.W[itime,:,j,i]
        Tprofile = self.T[itime,:,j,i]
        if overwrite:
            self.z = zprofile
            self.U = Uprofile
            self.V = Vprofile
            self.W = Wprofile
            self.T = Tprofile
        else:
            return dict(
                z=zprofile,
                U=Uprofile,
                V=Vprofile,
                W=Wprofile,
                T=Tprofile
            )

    def approx_z(self):
        self.zmean = self.z.mean(axis=(0,2,3))
        self.zstdev = self.z.std(axis=(0,2,3))
        return self.zmean

    def planar_average(self):
        """Horizontally average velocity and temperature fields

        Note: upwind fetch may skew the spatial average!
        """
        self.zmean = np.mean(self.z, axis=(0,2,3))
        self.Umean = np.mean(self.u, axis=(0,2,3))
        self.Vmean = np.mean(self.v, axis=(0,2,3))
        self.Wmean = np.mean(self.w, axis=(0,2,3))
        self.Tmean = np.mean(self.T, axis=(0,2,3))

