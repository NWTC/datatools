#
# Process WRF solution file
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import sys
import numpy as np

from netCDF4 import Dataset

g = 9.81

class WRFSolution(object):
    """Object to hold a single WRF solution snapshot"""

    def __init__(self,*args,**kwargs):
        verbose = kwargs.get('verbose',True)
        nclist = []
        Nfiles = len(args)
        self.filelist = args
        for i,fpath in enumerate(args):
            try:
                nc = Dataset(fpath)
            except OSError: # NetCDF: Unknown file format
                if verbose:
                    print('Skipped {:s} ({:d}/{:d})'.format(fpath,i+1,Nfiles))
            else:
                nclist.append(nc)
                if verbose:
                    print('Loaded {:s} ({:d}/{:d})'.format(fpath,i+1,Nfiles))

        if len(nclist) == 0:
            raise IOError('No WRF solution files found!')

        self.varlist = list(nc.variables)
        if verbose:
            print('  variables : ',self.varlist)

        # dimensions (unstaggered, i.e., face/cell-centered)
        self.Nt = nc.dimensions['time'].size
        self.Nx = nc.dimensions['west_east'].size 
        self.Ny = nc.dimensions['south_north'].size
        self.Nz = nc.dimensions['bottom_top'].size
        if verbose:
            print('  dimensions : ',self.Nt,self.Nx,self.Ny,self.Nz)

        # read variables and unstagger
        self._read_vars(nclist)


    def _read_vars(self,nclist):
        """Read all useful variables, shape==(NT,NZ,NY,NX)"""
        zlist = []
        Tlist = []
        Ulist = []
        Vlist = []
        Wlist = []
        for i,nc in enumerate(nclist):
            sys.stderr.write('\rReading vars from {:s}'.format(self.filelist[i]))
            Tlist.append( nc.variables['T'][:] )
            # staggered in x
            U = nc.variables['U'][:]
            Ulist.append( 0.5*(U[:,:,:,:-1] + U[:,:,:,1:]) )
            # staggered in y
            V = nc.variables['V'][:]
            Vlist.append( 0.5*(V[:,:,:-1,:] + V[:,:,1:,:]) )
            # staggered in z
            W = nc.variables['W'][:]
            PH = nc.variables['PH'][:]
            PHB = nc.variables['PHB'][:]
            Wlist.append( 0.5*(W[:,:-1,:,:] + W[:,1:,:,:]) )
            # calculate z == (ph + phb)/g
            zlist.append( 0.5*( PH[:,:-1,:,:] +  PH[:,1:,:,:] +
                               PHB[:,:-1,:,:] + PHB[:,1:,:,:] ) / g )

        print('Concatenating data arrays')
        self.U = np.concatenate(Ulist, axis=0)
        self.V = np.concatenate(Vlist, axis=0)
        self.W = np.concatenate(Wlist, axis=0)
        self.T = np.concatenate(Tlist, axis=0) + 300.0

        # calculate height AGL
        print('Calculating height AGL')
        self.z = np.concatenate(zlist, axis=0)
        if 'HGT' in self.varlist:
            # TODO: test this
            hgt = nclist[0].variables['HGT'][:]
            for i in range(self.Nx):
                for j in range(self.Ny):
                    self.z[:,i,j,:] -= hgt[i,j]

    
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
        if overwrite:
            self.z = self.z[itime,:,j,i]
            self.U = self.U[itime,:,j,i]
            self.V = self.V[itime,:,j,i]
            self.W = self.W[itime,:,j,i]
            self.T = self.T[itime,:,j,i]
        else:
            return {
                'z':self.z[itime,:,j,i],
                'U':self.U[itime,:,j,i],
                'V':self.V[itime,:,j,i],
                'W':self.W[itime,:,j,i],
                'T':self.T[itime,:,j,i]
            }

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

