
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

from ipywidgets import interactive, fixed #interact, interactive, fixed, interact_manual
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

    def plot(self, itime=-1, ax=None, speed_direction=False, **kwargs):
        if speed_direction:
            Uplot = np.sqrt(self.U[itime,:]**2 + self.V[itime,:]**2)
            Vplot = 180.0/np.pi * np.arctan2(-self.U[itime,:],-self.V[itime,:])
            Vplot[Vplot < 0] += 360.0
        else:
            Uplot = self.U[itime,:]
            Vplot = self.V[itime,:]
        if ax is None:
            fig,ax = plt.subplots(ncols=4,figsize=(10,4))
            fig.suptitle('t = {:.1f} s'.format(self.t[itime]))
        if len(self.z) > 1:
            ax[0].plot(Uplot, self.z, **kwargs)
            ax[1].plot(Vplot, self.z, **kwargs)
            ax[2].plot(self.W[itime,:], self.z, **kwargs)
        else:
            ax[0].axvline(Uplot, ls='--', **kwargs)
            ax[1].axvline(Vplot, ls='--', **kwargs)
            ax[2].axvline(self.W[itime,:], ls='--', **kwargs)
        if len(self.zT) > 1:
            ax[3].plot(self.T[itime,:], self.zT, **kwargs)
        else:
            ax[3].advline(self.T[itime,:], ls='--', **kwargs)
        if speed_direction:
            ax[0].set_xlabel(r'wind speed [m/s]')
            ax[1].set_xlabel(r'wind direction [deg]')
        else:
            ax[0].set_xlabel(r'$U$ [m/s]')
            ax[1].set_xlabel(r'$V$ [m/s]')
        ax[2].set_xlabel(r'$W$ [m/s]')
        ax[3].set_xlabel(r'$\theta$ [K]')
        ax[0].set_ylabel(r'$z$ [m]')

    def plot_over_time(self, **kwargs):
        """Plot profiles over time

        Keyword arguments
        -----------------
        speed_direction: bool, optional
            If True, plot wind speed and direction instead of U and V
            velocity components (default: False)
        convert_time: tuple, optional
            A single key/value pair where there key is the new time
            unit and the value is the time scaling factor
        time_range: tuple, optional
            Start/end times to plot (in converted time units); set to
            None for the beginning/end of the forcing data
        max_lines: int, optional
            Automatically calculate plotting interval if the number of
            snapshots to plot is greater than this number
        **kwargs: optional
            Passed to matplotlib.pyplot.plot
        """
        try:
            speeddir = kwargs.pop('speed_direction')
        except KeyError:
            speeddir = False

        try:
            tconv = kwargs.pop('convert_time')
        except KeyError:
            tconv = ('s', 1.0)
        tunit = tconv[0]
        tfac = tconv[1]

        try:
            trange = kwargs.pop('time_range')
        except KeyError:
            trange = (None,None)
        i0 = 0
        i1 = len(self.t)-1
        if trange[0] is not None:
            i0 = np.nonzero(self.t >= trange[0]/tfac)[0][0]
        if trange[1] is not None:
            i1 = np.nonzero(self.t <= trange[1]/tfac)[0][-1]
        subset = range(i0,i1+1)

        try:
            Nmax = kwargs.pop('max_lines')
        except KeyError:
            Nmax = 999999
        if len(subset) > Nmax:
            iskip = int(len(subset) / Nmax)
            subset = range(i0,i1+1,iskip)

        fig,ax = plt.subplots(ncols=4,figsize=(10,4))
        colors = cm.get_cmap(series_colormap)
        frac = (subset - i0) / (i1 - i0)
        for i, itime in enumerate(subset):
            ti = self.t[itime]
            col = colors(frac[i])
            label = ''
            if (i == 0) or (i == len(subset)-1):
                label = '{:.1f} {:s}'.format(ti*tfac,tunit)
            kwargs['color'] = col
            kwargs['label'] = label
            self.plot(itime=itime, ax=ax, speed_direction=speeddir, **kwargs)
        ax[-1].legend(loc='upper left',bbox_to_anchor=(1.05,1.0))

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

    def _calculate_Tgrad_upper(self):
        iz = np.nonzero(self.z >= self.inversion_top.value)[0][0]
        dz = self.z[-1] - self.z[iz]
        assert(dz > 0)
        return (self.T[:,-1] - self.T[:,iz]) / dz

    def editor_plot(self,**kwargs):
        edits = self.editor.kwargs
        plot_kwargs = kwargs.get('plot_kwargs',{})
        if 'max_lines' not in plot_kwargs.keys():
            plot_kwargs['max_lines'] = int((self.t[-1]-self.t[0])/3600.0) + 1

        conv = 3600.0
        TGradUpper = None

        # update controls
        if edits['enforce_lapse_rate']:
            self.lapse_rate.disabled = False
            self.inversion_top.disabled = False
            TGradUpper = self._calculate_Tgrad_upper()
            self.lapse_rate.min = np.min(TGradUpper)
            self.lapse_rate.max = np.max(TGradUpper)
        else:
            self.lapse_rate.disabled = True
            self.inversion_top.disabled = True

        # first time setup
        if not hasattr(self, '_Uext'):
            _,Nz = self.U.shape
            _,NzT = self.T.shape
            self._text = np.zeros((self.Nt+2))
            self._Uext = np.zeros((self.Nt+2,Nz))
            self._Vext = np.zeros((self.Nt+2,Nz))
            self._Wext = np.zeros((self.Nt+2,Nz))
            self._Text = np.zeros((self.Nt+2,NzT))
            self._lapse_correction = np.zeros((self.Nt+2,NzT))
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

        # extrapolate temperature profile
        if edits['enforce_lapse_rate']:
            iz = np.nonzero(self.z >= self.inversion_top.value)[0][0]
            self._lapse_correction[:,:] = 0.0
            delta = self.z[iz:] - self.z[iz]
            for it in range(len(self._text)):
                self._lapse_correction[it,iz:] = -self._Text[it,iz:] \
                                                + self._Text[it,iz] \
                                                + delta*self.lapse_rate.value

        # plot extended values
        self.t = self._text
        self.U = self._Uext
        self.V = self._Vext
        self.W = self._Wext
        if edits['enforce_lapse_rate']:
            self.T = self._Text + self._lapse_correction
        else:
            self.T = self._Text

        self.plot_over_time(**plot_kwargs)

        if edits['enforce_lapse_rate']:
            plt.gcf().get_axes()[3].axhline(edits['inversion_top'],
                                            ls='--',color='k')

        # restore actual values
        self.t = self._tsave
        self.U = self._Usave
        self.V = self._Vsave
        self.W = self._Wsave
        self.T = self._Tsave

    def edit(self,**kwargs):
        """Basic controls for manipulating the source terms"""
        self.start_hrs = widgets.BoundedFloatText(value=0.0,min=-999,max=0.0,step=0.25,
                                                  description='hours')
        self.end_hrs = widgets.BoundedFloatText(value=0.0,min=0.0,max=999,step=0.25,
                                                description='hours')

        self.inversion_top = widgets.FloatSlider(min=self.z[1], max=self.z[-2],
                                                 value=self.z[-2],
                                                 step=self.z[2]-self.z[0],
                                                 readout_format='.1f',
                                                 disabled=True)
        TGradUpper = self._calculate_Tgrad_upper()
        self.lapse_rate = widgets.FloatSlider(min=np.min(TGradUpper),
                                              max=np.max(TGradUpper),
                                              step=0.0001,
                                              value=np.mean(TGradUpper),
                                              readout_format='.3g',
                                              disabled=True)
        self.editor = interactive(self.editor_plot,
                                  plot_kwargs=fixed(kwargs),
                                  mom_start=['extend constant','extrapolate'],
                                  temp_start=['extend constant','extrapolate'],
                                  start_hrs=self.start_hrs,
                                  mom_end=['extend constant','extrapolate'],
                                  temp_end=['extend constant','extrapolate'],
                                  end_hrs=self.end_hrs,
                                  enforce_lapse_rate=False,
                                  lapse_rate=self.lapse_rate,
                                  inversion_top=self.inversion_top
                                  )
        display(self.editor)

    def reset(self):
        delattr(self, '_Uext') # flag for update
        delattr(self, '_Vext')
        delattr(self, '_Wext')
        delattr(self, '_Text')
        delattr(self, '_text')
        delattr(self, '_lapse_correction')
        delattr(self, '_tsave')
        delattr(self, '_Usave')
        delattr(self, '_Vsave')
        delattr(self, '_Wsave')
        delattr(self, '_Tsave')

    def save_edits(self):
        edits = self.editor.kwargs
        istart, iend = None, None
        if edits['start_hrs'] == 0:
            istart = 1
        if edits['end_hrs'] == 0:
            iend = -1
        inrange = slice(istart,iend)
        self.U = self._Uext[inrange,:]
        self.V = self._Vext[inrange,:]
        self.W = self._Wext[inrange,:]
        self.T = self._Text[inrange,:]
        if edits['enforce_lapse_rate']:
            self.T += self._lapse_correction[inrange,:]
        self.t = self._text[inrange]
        toffset = -self.t[0]
        print('shifting time by {:.1f} s'.format(toffset))
        self.t += toffset
        self.Nt = len(self.t)
        self.reset()

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

