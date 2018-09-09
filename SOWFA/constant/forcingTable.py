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

from datatools.openfoam_util import read_all_defs, of_list, of_listlist

from ipywidgets import interactive #interact, interactive, fixed, interact_manual
import ipywidgets as widgets
from IPython.display import display

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
        self.z = np.array(heights)
        self.zT = self.z
        self.t = np.array(times)
        if any([ profile is not None for profile in [U,V,W,T]]):
            # if specified this way, assume all profiles have the same heights
            for profile in [U,V,W,T]:
                if profile is not None:
                    break
            Nt, Nz = profile.shape
            self.Nt = Nt
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
        Nt = self.Nt
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
            self.t = np.array(t)
        if self.Nt == 1:
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
        try:
            tconv = kwargs.pop('convert_time')
        except KeyError:
            tconv = {'s': 1.0}
        tunit = list(tconv.keys())[0]
        tfac = tconv[tunit]

        fig,ax = plt.subplots(ncols=4,figsize=(10,4))
        colors = cm.get_cmap(series_colormap)
        frac = (self.t - self.t[0]) / (self.t[-1] - self.t[0])
        for itime, ti in enumerate(self.t):
            col = colors(frac[itime])
            label = ''
            if (itime == 0) or (itime == len(self.t)-1):
                label = '{:.1f} {:s}'.format(ti*tfac,tunit)
            kwargs['color'] = col
            kwargs['label'] = label
            self.plot(itime=itime, ax=ax, **kwargs)
        ax[0].legend(loc='best')

    def extrapolate(self,field,time,order=1):
        from scipy.interpolate import InterpolatedUnivariateSpline
        U = getattr(self,field)
        _,Nz = U.shape
        interpname = '_' + field + 'interp'
        try:
            Ufn = getattr(self, interpname)
        except AttributeError:
            Ufn = []
            for iz in range(Nz):
                Ufn.append(InterpolatedUnivariateSpline(self.t, U[:,iz], k=order))
            setattr(self, interpname, Ufn)
        Uextrap = np.zeros((Nz))
        for iz in range(Nz):
            Uextrap[iz] = Ufn[iz](time)
        return Uextrap

    def editor_plot(self,**kwargs):
        edits = self.editor.kwargs
        conv = 3600.0
        if not hasattr(self, '_Uext'):
            _,Nz = self.U.shape
            _,NzT = self.T.shape
            self._text = np.zeros((self.Nt+2))
            self._Uext = np.zeros((self.Nt+2,Nz))
            self._Vext = np.zeros((self.Nt+2,Nz))
            self._Wext = np.zeros((self.Nt+2,Nz))
            self._Text = np.zeros((self.Nt+2,NzT))
            self._text[1:-1] = self.t
            self._Uext[1:-1,:] = self.U
            self._Vext[1:-1,:] = self.V
            self._Wext[1:-1,:] = self.W
            self._Text[1:-1,:] = self.T
            self._tsave = self.t.copy()
            self._Usave = self.U.copy()
            self._Vsave = self.V.copy()
            self._Wsave = self.W.copy()
            self._Tsave = self.T.copy()
        # get updated times
        t0 = self.t[0] + edits['start_hrs'] * conv
        t1 = self.t[-1] + edits['end_hrs'] * conv
        self._text[0] = t0
        self._text[-1] = t1
        # extrapolate
        if edits['mom_start'] == 'extrapolate':
            self._Uext[0,:] = self.extrapolate('U',t0)
            self._Vext[0,:] = self.extrapolate('V',t0)
            self._Wext[0,:] = self.extrapolate('W',t0)
        else:
            self._Uext[0,:] = self._Usave[0,:]
            self._Vext[0,:] = self._Vsave[0,:]
            self._Wext[0,:] = self._Wsave[0,:]
        if edits['temp_start'] == 'extrapolate':
            self._Text[0,:] = self.extrapolate('T',t0)
        else:
            self._Text[0,:] = self._Tsave[0,:]
        if edits['mom_end'] == 'extrapolate':
            self._Uext[-1,:] = self.extrapolate('U',t1)
            self._Vext[-1,:] = self.extrapolate('V',t1)
            self._Wext[-1,:] = self.extrapolate('W',t1)
        else:
            self._Uext[-1,:] = self._Usave[-1,:]
            self._Vext[-1,:] = self._Vsave[-1,:]
            self._Wext[-1,:] = self._Wsave[-1,:]
        if edits['temp_end'] == 'extrapolate':
            self._Text[-1,:] = self.extrapolate('T',t1)
        else:
            self._Text[-1,:] = self._Tsave[-1,:]
        # plot extrapolated values
        self.t = self._text
        self.U = self._Uext
        self.V = self._Vext
        self.W = self._Wext
        self.T = self._Text
        self.plot_over_time(convert_time={'h':1.0/conv})
        # restore actual values
        self.t = self._tsave
        self.U = self._Usave
        self.V = self._Vsave
        self.W = self._Wsave
        self.T = self._Tsave

    def edit(self):
        """Basic controls for manipulating the source terms"""
        self.start_hrs = widgets.BoundedFloatText(value=0.0,min=-999,max=0.0,step=0.25,
                                                  description='hours')
        self.end_hrs = widgets.BoundedFloatText(value=0.0,min=0.0,max=999,step=0.25,
                                                description='hours')
        self.editor = interactive(self.editor_plot,
                                  mom_start=['extend constant','extrapolate'],
                                  temp_start=['extend constant','extrapolate'],
                                  start_hrs=self.start_hrs,
                                  mom_end=['extend constant','extrapolate'],
                                  temp_end=['extend constant','extrapolate'],
                                  end_hrs=self.end_hrs)
        display(self.editor)

    def save_edits(self):
        edits = self.editor.kwargs
        istart, iend = None, None
        if edits['start_hrs'] == 0:
            istart = 1
        if edits['end_hrs'] == 0:
            iend = -1
        inrange = slice(istart,iend)
        self.U = self._Uext[inrange]
        self.V = self._Vext[inrange]
        self.W = self._Wext[inrange]
        self.T = self._Text[inrange]
        self.t = self._text[inrange]
        toffset = -self.t[0]
        print('shifting time by {:.1f} s'.format(toffset))
        self.t += toffset
        self.Nt = len(self.t)
        delattr(self, '_Uext') # flag for update

    def to_csv(self,fname):
        alltimesM = [ ti for ti in self.t for _ in range(len(self.z)) ]
        df = pd.DataFrame(index=alltimesM)
        df['z'] = np.tile(self.z, self.Nt)
        df['U'] = self.U.ravel()
        df['V'] = self.V.ravel()
        df['W'] = self.W.ravel()
        if not self.separate_heights:
            df['T'] = self.T.ravel()
        else:
            alltimesT = [ ti for ti in self.t for _ in range(len(self.zT)) ]
            dfT = pd.DataFrame(index=alltimesT)
            dfT['z'] = np.tile(self.zT, self.Nt)
            dfT['T'] = self.T.ravel()
            df = pd.concat([df,dfT],sort=False) # suppress warning, don't sort columns
        df.to_csv(fname)
        print('Wrote '+fname)

    def read_csv(self,fname):
        """Read forcingTable.csv generated by this module"""
        df = pd.read_csv(fname,index_col=0)
        self.t = df.index.unique()
        self.z = df.loc[~pd.isna(df['U']),'z'].unique()
        self.zT = df.loc[~pd.isna(df['T']),'z'].unique()
        self.Nt = len(self.t)
        Nz = len(self.z)
        NzT = len(self.zT)
        U = df.loc[~pd.isna(df['U']),'U']
        V = df.loc[~pd.isna(df['V']),'V']
        W = df.loc[~pd.isna(df['W']),'W']
        T = df.loc[~pd.isna(df['T']),'T']
        self.U = U.values.reshape((self.Nt,Nz))
        self.V = V.values.reshape((self.Nt,Nz))
        self.W = W.values.reshape((self.Nt,Nz))
        self.T = T.values.reshape((self.Nt,NzT))
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
            newdata = read_all_defs(fname,verbose=False)
            if data is None:
                data = newdata
            else:
                data.update(newdata)

        self.t = data[U][:,0]
        assert(np.all(self.t==data[V][:,0]))
        assert(np.all(self.t==data[W][:,0]))
        assert(np.all(self.t==data[T][:,0]))
        self.Nt = len(self.t)

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
        
    def write(self,fname='forcingTable',
            z='sourceHeightsMomentum',
            U='sourceTableMomentumX',
            V='sourceTableMomentumY',
            W='sourceTableMomentumZ',
            zT='sourceHeightsTemperature',
            T='sourceTableTemperature',
            ):
        """Write out forcingTable to be included in a SOWFA simulation
        (e.g., within constant/ABLProperties). The specified variable
        names are the data arrays expected by ABLSolver
        """
        # prepend the time for each row of source terms
        timeU = np.concatenate((self.t[:,np.newaxis], self.U), axis=1)
        timeV = np.concatenate((self.t[:,np.newaxis], self.V), axis=1)
        timeW = np.concatenate((self.t[:,np.newaxis], self.W), axis=1)
        timeT = np.concatenate((self.t[:,np.newaxis], self.T), axis=1)
        # now dump everything out
        with open(fname,'w') as f:
            f.write(of_list(z,self.z))
            f.write(of_listlist(U,timeU))
            f.write(of_listlist(V,timeV))
            f.write(of_listlist(W,timeW))
            f.write(of_list(zT,self.zT))
            f.write(of_listlist(T,timeT))
        print('Wrote '+fname)

