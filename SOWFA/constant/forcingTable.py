#
# Module for read/writing SOWFA source terms, e.g.:
#   constant/forcingTable
# stored in OpenFOAM's table format, referenced by:
#   constant/ABLProperties.
#
# This is also used by the WRF.viz module for WRF-to-SOWFA
# mesoscale-to-microscale coupling.
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd

from datatools import openfoam_util

series_colormap = 'viridis'

class ForcingTable(object):
    def __init__(self,heights=None,times=None,U=None,V=None,W=None,T=None):
        """Create ForcingTable object. Velocity (U,V,W) and potential
        temperature (T) time-height profiles may be specified at this time
        or read in later
        """
        self.U = U
        self.V = V
        self.W = W
        self.T = T
        self.separate_heights = False
        self.z = heights
        self.zT = self.z
        self.t = times
        if any([ profile is not None for profile in [U,V,W,T]]):
            # if specified this way, assume all profiles have the same heights
            for profile in [U,V,W,T]:
                if profile is not None:
                    break
            Nt, Nz = profile.shape
            assert(Nt == len(self.t))
            if U is None:
                self.U = np.zeros((Nt,Nz))
            if V is None:
                self.V = np.zeros((Nt,Nz))
            if W is None:
                self.W = np.zeros((Nt,Nz))
            if T is None:
                self.T = np.zeros((Nt,Nz))
            self._validate()

    def regularize_heights(self, z):
        """Given heights from WRF, heights are not guarnateed to be
        constant. This removes the time dimension from height and
        interpolates all field variables to the specified z
        """
        Nz = len(z)
        Nt = len(self.t)
        Unew = np.zeros((Nt,Nz))
        Vnew = np.zeros((Nt,Nz))
        Wnew = np.zeros((Nt,Nz))
        Tnew = np.zeros((Nt,Nz))
        for itime,ti in enumerate(self.t):
            if len(self.z) > 1:
                print('interpolating time {:d}: t={:s}'.format(itime,str(ti)))    
                Unew[itime,:] = np.interp(z, self.z[itime,:], self.U[itime,:])
                Vnew[itime,:] = np.interp(z, self.z[itime,:], self.V[itime,:])
                Wnew[itime,:] = np.interp(z, self.z[itime,:], self.W[itime,:])
            else:
                Unew[itime,:] = self.U[itime,0]
                Vnew[itime,:] = self.V[itime,0]
                Wnew[itime,:] = self.W[itime,0]
            if len(self.zT) > 1:
                if self.separate_heights:
                    Tnew[itime,:] = np.interp(z, self.zT[itime,:], self.T[itime,:])
                else:
                    Tnew[itime,:] = np.interp(z, self.z[itime,:], self.T[itime,:])
            else:
                Tnew[itime,:] = self.T[itime,0]
        self.z = z
        self.zT = z
        self.U = Unew
        self.V = Vnew
        self.W = Wnew
        self.T = Tnew
        
    def _validate(self):
        assert(self.z is not None)
        assert(self.U.shape == self.V.shape == self.W.shape)
        if isinstance(self.t[0], datetime):
            # assume all times are datetime (or pandas timestamp) objects
            print('converting time stamps to seconds')
            t = []
            for ti in self.t:
                dt = ti - self.t[0]
                t.append(dt.total_seconds())
            self.t = t
        if len(self.t) == 1:
            # duplicate profile for t=TLARGE so that we have constant source terms
            print('duplicating time 0 for constant profile')
            self.t = [self.t[0], self.t[0] + 999999.0]
            self.U = np.tile(self.U,[2,1])
            self.V = np.tile(self.V,[2,1])
            self.W = np.tile(self.W,[2,1])
            self.T = np.tile(self.T,[2,1])

    def plot(self, itime=-1, ax=None, **kwargs):
        if ax is None:
            fig,ax = plt.subplots(ncols=4,figsize=(10,4))
            fig.suptitle('t = {:.1f} s'.format(self.t[itime]))
        if len(self.z) > 1:
            ax[0].plot(self.U[itime,:], self.z, **kwargs)
            ax[1].plot(self.V[itime,:], self.z, **kwargs)
            ax[2].plot(self.W[itime,:], self.z, **kwargs)
        else:
            ax[0].axvline(self.U[itime,:], ls='--', **kwargs)
            ax[1].axvline(self.V[itime,:], ls='--', **kwargs)
            ax[2].axvline(self.W[itime,:], ls='--', **kwargs)
        if len(self.zT) > 1:
            ax[3].plot(self.T[itime,:], self.zT, **kwargs)
        else:
            ax[3].advline(self.T[itime,:], ls='--', **kwargs)
        ax[0].set_xlabel(r'$U$ [m/s]')
        ax[1].set_xlabel(r'$V$ [m/s]')
        ax[2].set_xlabel(r'$W$ [m/s]')
        ax[3].set_xlabel(r'$\theta$ [K]')
        ax[0].set_ylabel(r'$z$ [m]')

    def plot_over_time(self, **kwargs):
        fig,ax = plt.subplots(ncols=4,figsize=(10,4))
        colors = cm.get_cmap(series_colormap)
        for itime, ti in enumerate(self.t):
            col = colors(float(itime)/(len(self.t)-1))
            label = ''
            if (itime == 0) or (itime == len(self.t)-1):
                label = 't = {:.1f} s'.format(ti)
            kwargs['color'] = col
            kwargs['label'] = label
            self.plot(itime=itime, ax=ax, **kwargs)
        ax[0].legend(loc='best')

    def to_csv(self,fname):
        alltimesM = [ ti for ti in self.t for _ in range(len(self.z)) ]
        df = pd.DataFrame(index=alltimesM)
        df['z'] = np.tile(self.z, len(self.t))
        df['U'] = self.U.ravel()
        df['V'] = self.V.ravel()
        df['W'] = self.W.ravel()
        if not self.separate_heights:
            df['T'] = self.T.ravel()
        else:
            alltimesT = [ ti for ti in self.t for _ in range(len(self.zT)) ]
            dfT = pd.DataFrame(index=alltimesT)
            dfT['z'] = np.tile(self.zT, len(self.t))
            dfT['T'] = self.T.ravel()
            df = pd.concat([df,dfT],sort=False) # suppress warning, don't sort columns
        df.to_csv(fname)

    def read_csv(self,fname):
        """Read forcingTable.csv generated by this module"""
        df = pd.read_csv(fname,index_col=0)
        self.t = df.index.unique()
        self.z = df.loc[~pd.isna(df['U']),'z'].unique()
        self.zT = df.loc[~pd.isna(df['T']),'z'].unique()
        Nt = len(self.t)
        Nz = len(self.z)
        NzT = len(self.zT)
        U = df.loc[~pd.isna(df['U']),'U']
        V = df.loc[~pd.isna(df['V']),'V']
        W = df.loc[~pd.isna(df['W']),'W']
        T = df.loc[~pd.isna(df['T']),'T']
        self.U = U.values.reshape((Nt,Nz))
        self.V = V.values.reshape((Nt,Nz))
        self.W = W.values.reshape((Nt,Nz))
        self.T = T.values.reshape((Nt,NzT))
        self._validate()

    def read_openfoam_ascii(self,*args,
            z='sourceHeightsMomentum',
            U='sourceTableMomentumX',
            V='sourceTableMomentumY',
            W='sourceTableMomentumZ',
            zT='sourceHeightsTemperature',
            T='sourceTableTemperature',
            ):
        """Read forcingTable(s) suitable to be included in a SOWFA 
        simulation (e.g., within constant/ABLProperties). The specified
        variable names are read from OpenFOAM tables.
        """
        data = None    
        for fname in args:
            newdata = openfoam_util.read_all_tables(fname)
            if data is None:
                data = newdata
            else:
                data.update(newdata)

        self.t = data[U][:,0]
        assert(np.all(self.t==data[V][:,0]))
        assert(np.all(self.t==data[W][:,0]))
        assert(np.all(self.t==data[T][:,0]))

        self.z = data[z]
        self.U = data[U][:,1:]
        self.V = data[V][:,1:]
        self.W = data[W][:,1:]
        self.zT = data[zT]
        self.T = data[T][:,1:]
        if (not len(self.z) == len(self.zT)) \
                or ~np.all(self.z == self.zT):
            self.separate_heights = True

        self._validate()
        
    #def to_openfoam(self,fname):

