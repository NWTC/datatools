#
# WRF visualization helper module
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import os
from ipywidgets import interactive #interact, interactive, fixed, interact_manual
import ipywidgets as widgets

import xarray
import numpy as np
import matplotlib.pyplot as plt

g = 9.81

def _reorient(M):
    #return np.fliplr(M).T
    return np.fliplr(M)

class Visualization2D(object):

    def __init__(self,*args,**kwargs):
        tdim = kwargs.get('tdim','time') # time / file aggregation dimension ("aggdim")
        xdim = kwargs.get('xdim','west_east') # unstaggered by default
        ydim = kwargs.get('ydim','south_north') # unstaggered by default
        zdim = kwargs.get('zdim','bottom_top') # unstaggered by default
        plane = kwargs.get('plane','z') # 2D plane normal direction

        """Load a series of netcdf files provided by args"""
        if len(args) > 0:
            inputfiles = args
        else:
            inputfiles = os.listdir('.')
        filelist = []
        for fpath in inputfiles:
            if not os.path.isfile(fpath):
                continue
            try:
                xarray.open_dataset(fpath)
            except (IOError,OSError):
                # NetCDF: Unknown file format
                continue
            else:
                filelist.append(fpath)
        filelist.sort()
        self.filelist = filelist
        self.data = xarray.open_mfdataset(filelist, concat_dim=tdim)

        """Set up dimensions"""
        self.plane = plane
        assert(plane in ['x','y','z'])
        self.Ntimes = self.data.dims[tdim]
        self.Nx = self.data.dims[xdim]
        self.Ny = self.data.dims[ydim]
        self.Nz = self.data.dims[zdim]
        if plane == 'x':
            self.N = self.Nx
        elif plane == 'y':
            self.N = self.Ny
        elif plane == 'z':
            self.N = self.Nz

        """Set up field variables"""
        # unstaggered
        self.T = self.data.variables['T'][:] + 300.0
        # staggered in x
        U = self.data.variables['U'][:]
        self.U = 0.5*(U[:,:,:,:-1] + U[:,:,:,1:])
        # staggered in y
        V = self.data.variables['V'][:]
        self.V = 0.5*(V[:,:,:-1,:] + V[:,:,1:,:])
        # staggered in z
        W = self.data.variables['W'][:]
        PH = self.data.variables['PH'][:]
        PHB = self.data.variables['PHB'][:]
        self.W = 0.5*(W[:,:-1,:,:] + W[:,1:,:,:])
        # calculate z = (ph + phb)/g
        self.z = 0.5*( PH[:,:-1,:,:] +  PH[:,1:,:,:] +
                      PHB[:,:-1,:,:] + PHB[:,1:,:,:] ) / g
        # other variables
        #self.Umag = np.sqrt(self.U**2 + self.V**2 + self.W**2) # can cause a memory error
        
    def __repr__(self):
        s = str(self.Ntimes) + ' times read:\n'
        for fpath in self.filelist[:3]:
            s += '  ' + os.path.split(fpath)[-1] + '\n'
        if len(self.filelist) > 5:
            s += '  ...\n'
        for fpath in self.filelist[-2:]:
            s += '  ' + os.path.split(fpath)[-1] + '\n'
        s += 'Dimensions: ({:d}, {:d}, {:d})'.format(self.Nx,
                                                     self.Ny,
                                                     self.Nz)
        return s

    def plot(self,field='Umag',time=0,index=0):
        """Jupyter notebook usage:
            iplot = interactive(viz.plot,
                                plane=['z','y','x'],
                                field=['U','V','W','T'],
                                time=(0,viz.Ntimes-1),
                                index=(0,viz.N-1))
            iplot
        """
        assert((time >= 0) and (time < self.Ntimes))
        plt.figure(1)
        if self.plane == 'z':
            assert((index >= 0) and (index < self.Nz))
            U = getattr(self,field)
            U = U[time,index,:,:]
            cont = plt.imshow(_reorient(U))
        else:
            print(self.plane,'not supported')
        cbar = plt.colorbar(cont)
        cbar.set_label(field)
        plt.show()
        

