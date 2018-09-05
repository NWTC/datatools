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
import pandas as pd

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
        self.z = heights
        self.t = times
        if any([ profile is not None for profile in [U,V,W,T]]):
            for profile in [U,V,W,T]:
                if profile is not None:
                    break
            self.Nt, self.Nz = profile.shape
            assert(self.Nt == len(self.t))
            if U is None:
                self.U = np.zeros((self.Nt, self.Nz))
            if V is None:
                self.V = np.zeros((self.Nt, self.Nz))
            if W is None:
                self.W = np.zeros((self.Nt, self.Nz))
            if T is None:
                self.T = np.zeros((self.Nt, self.Nz))
            self._validate()

    def regularize_heights(self, z):
        """Given heights from WRF, heights are not guarnateed to be
        constant. This removes the time dimension from height and
        interpolates all field variables to the specified z
        """
        self.Nz = len(z)
        Unew = np.zeros((self.Nt,self.Nz))
        Vnew = np.zeros((self.Nt,self.Nz))
        Wnew = np.zeros((self.Nt,self.Nz))
        Tnew = np.zeros((self.Nt,self.Nz))
        for itime,ti in enumerate(self.t):
            print('interpolating time {:d}: t={:s}'.format(itime,str(ti)))    
            Unew[itime,:] = np.interp(z, self.z[itime,:], self.U[itime,:])
            Vnew[itime,:] = np.interp(z, self.z[itime,:], self.V[itime,:])
            Wnew[itime,:] = np.interp(z, self.z[itime,:], self.W[itime,:])
            Tnew[itime,:] = np.interp(z, self.z[itime,:], self.T[itime,:])
        self.z = z
        self.U = Unew
        self.V = Vnew
        self.W = Wnew
        self.T = Tnew
        
    def _validate(self):
        assert(self.z is not None)
        assert(self.U.shape == self.V.shape == self.W.shape == self.T.shape)
        assert(self.U.shape == (self.Nt, self.Nz))
        if isinstance(self.t[0], datetime):
            # assume all times are datetime (or pandas timestamp) objects
            print('converting time stamps to seconds')
            t = []
            for ti in self.t:
                dt = ti - self.t[0]
                t.append(dt.total_seconds())
            self.t = t
        if self.Nt == 1:
            # duplicate profile for t=TLARGE so that we have constant source terms
            print('duplicating time 0 for constant profile')
            self.Nt = 2
            self.t = [self.t[0], self.t[0] + 999999.0]
            self.U = np.tile(self.U,[2,1])
            self.V = np.tile(self.V,[2,1])
            self.W = np.tile(self.W,[2,1])
            self.T = np.tile(self.T,[2,1])

    def to_csv(self,fname):
        alltimes = [ ti for ti in self.t for _ in range(self.Nz) ]
        df = pd.DataFrame(index=alltimes)
        df['z'] = np.tile(self.z, self.Nt)
        df['U'] = self.U.ravel()
        df['V'] = self.V.ravel()
        df['W'] = self.W.ravel()
        df['T'] = self.T.ravel()
        df.to_csv(fname)

    def read_csv(self,fname):
        df = pd.read_csv(fname,index_col=0)
        self.t = df.index.unique()
        self.z = df['z'].unique()
        self.Nt = len(self.t)
        self.Nz = len(self.z)
        self.U = df['U'].values.reshape((self.Nt, self.Nz))
        self.V = df['V'].values.reshape((self.Nt, self.Nz))
        self.W = df['W'].values.reshape((self.Nt, self.Nz))
        self.T = df['T'].values.reshape((self.Nt, self.Nz))
        
    #def to_openfoam(self,fname):

