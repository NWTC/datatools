"""
For processing SOWFA cell averages in postProcessing/averaging
Based on original SOWFA-Tools in MATLAB

written by Eliot Quon (eliot.quon@nrel.gov)

Sample usage:

    from SOWFA.postProcessing.averaging import read

    # read all time directories in current working directory
    averagingData = read()

    # read '0' and '1000' in current working directory
    averagingData = read( 0, 1000 )

    # read all time directories in specified directory
    averagingData = read('caseX/postProcessing/averaging')

    # read specified time directories
    averagingData = read('caseX/postProcessing/averaging/0',
                        'caseX/postProcessing/averaging/1000')

"""
from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

all_vars = [
        'U_mean','V_mean','W_mean','T_mean',
        'uu_mean', 'vv_mean', 'ww_mean', 'uv_mean', 'uw_mean', 'vw_mean',
        'Tw_mean',
        'R11_mean','R22_mean','R33_mean','R12_mean','R13_mean','R23_mean',
        'q1_mean','q2_mean','q3_mean',
        ]

class PlanarAverages(object):

    def __init__(self,*args,**kwargs):
        """Find and process all time directories"""
        self._processed = []
        self.hLevelsCell = None
        self.simTimeDirs = [] # output time names
        self.simStartTimes = [] # start or restart simulation times
        self.imax = None # for truncating time series
        self.rotated = False

        if len(args)==0:
            args = os.listdir('.')
            if 'averaging' in args: args = ['averaging']

        # find results
        for arg in args:
            try:
                if not os.path.isdir(arg): continue
            except TypeError: # a number was specified
                arg = '{:g}'.format(arg)
                if not os.path.isdir(arg): continue

            if arg[-1] == os.sep: arg = arg[:-1] # strip trailing slash

            listing = os.listdir(arg)
            if 'hLevelsCell' in listing:
                # an output (time) directory was directly specified
                self.simTimeDirs.append(arg)
                try:
                    timeDirName = os.path.split(arg)[1] # final part of path after last slash
                    dirTime = float(timeDirName)

                except: # specified results dir is not a number
                    dirTime = -1
                self.simStartTimes.append(dirTime)
            elif not arg.startswith('boundaryData'):
                print('Checking directory:',arg) #,'with',listing
                # specified a directory containing output (time) subdirectories
                for dirname in listing:
                    if not os.path.isdir(arg+os.sep+dirname): continue
                    #print('  checking subdirectory',dirname)
                    try:
                        startTime = float(dirname)
                        if 'hLevelsCell' in os.listdir(arg+os.sep+dirname):
                            self.simTimeDirs.append( arg+os.sep+dirname )
                            self.simStartTimes.append( startTime )
                    except ValueError:
                        # dirname is not a number
                        pass

        # sort results
        self.simTimeDirs = [ x[1] for x in sorted(zip(self.simStartTimes,self.simTimeDirs)) ]
        self.simStartTimes.sort()

        print('Simulation (re)start times:',self.simStartTimes)

        # process all output dirs
        #for idir,tdir in enumerate( self.simTimeDirs ):
        #    print('Processing',tdir)
        #    self._process(tdir)
        if len(self.simTimeDirs) > 0:
            self._processdirs( self.simTimeDirs, **kwargs )
        else:
            print('No averaging time directories found!')
    
        self._trim_series_if_needed()


    def __repr__(self):
        s = 'SOWFA postProcessing: averaging data'
        for t,d in zip(self.simStartTimes,self.simTimeDirs):
            fullpath = os.path.realpath(os.curdir) + os.sep + d
            s += '\n  {:f}\t{:s}'.format(t,fullpath)
        s += '\ntimes: {:d} [{:f},{:f}]'.format(len(self.t),self.t[0],self.t[-1])
        s += '\nheights: {:d} [{:f},{:f}]'.format(len(self.hLevelsCell),self.hLevelsCell[0],self.hLevelsCell[-1])
        return s

    def get_vars_if_needed(self,*args,**kwargs):
        """Read in specified list of variables"""
        varList = []
        reread = kwargs.pop('reread',False)
        if reread:
            for var in args:
                if var in self._processed: self._processed.remove(var)
                varList.append(var)
            print('Rereading variables',varList)
        else:
            for var in args:
                if var not in self._processed: varList.append(var)
        if len(varList)==0: return

        kwargs['varList'] = varList
        self._processdirs( self.simTimeDirs, **kwargs )

    def _processdirs(self,
                     tdirList,
                     varList=['U_mean','V_mean','W_mean','T_mean'],
                     trimOverlap=True
                    ):
        """Reads all files within an averaging output time directory,
        presumably containing hLevelsCell and other cell-averaged
        quantities. An object attribute corresponding to the averaged
        output name is updated, e.g.:
            ${timeDir}/U_mean is appended to the array self.U_mean

        Typically, objects have shape (Nt,Nz).
        """
        outputs = []
        if isinstance( varList, (str,) ):
            if varList.lower()=='all':
                # special case: read all vars
                allOutputs = os.listdir(tdirList[0])
                for field in allOutputs:
                    if field=='hLevelsCell':
                        continue
                    else:
                        outputs.append( field )
            else: # specified single var
                outputs = [varList]
        else: # specified list
            outputs = varList

        # process hLevelsCell first, verify we have the same cells
        with open(tdirList[0]+os.sep+'hLevelsCell','r') as f:
            line = f.readline()
        self.hLevelsCell = np.array([ float(val) for val in line.split() ])

        # check that we have the same amount of data in all fields
        for tdir in tdirList:
            Nlines = []
            for field in outputs:
                output = tdir + os.sep + field
                if not os.path.isfile(output):
                    print('Error:',output,'not found')
                    return

                with open(output,'r') as f:
                    for i,line in enumerate(f): pass
                    Nlines.append(i+1)
                    line = line.split() # check final line for the right number of values
                    if not len(line) == len(self.hLevelsCell)+2: # t,dt,f_1,f_2,...,f_N for N heights
                        print('z',z)
                        print('line',line)
                        print('Error: number of output points inconsistent with hLevelsCell in', output)
                        return
            if not np.min(Nlines) == np.max(Nlines):
                print('Warning: number of output times do not match in all files')
        N = Nlines[0]

        # NOW process all data
        selected = []
        for field in outputs:
            arrays = [ np.loadtxt( tdir + os.sep + field ) for tdir in tdirList ]

            # combine into a single array and trim end of time series
            # (because simulations that are still running can have different
            # array lengths)
            newdata = np.concatenate(arrays)[:self.imax,:]

            # get rid of overlapped data for restarts
            if trimOverlap:
                if len(selected) == 0:
                    # create array mask
                    tpart = [ array[:,0] for array in arrays ]
                    for ipart,tcutoff in enumerate(self.simStartTimes[1:]):
                        selectedpart = np.ones(len(tpart[ipart]),dtype=bool)
                        try:
                            iend = np.nonzero(tpart[ipart] > tcutoff)[0][0]
                        except IndexError:
                            # clean restart
                            pass
                        else:
                            # previous simulation didn't finish; overlapped data
                            selectedpart[iend:] = False 
                        selected.append(selectedpart)
                    # last / currently running part
                    selected.append(np.ones(len(tpart[-1]),dtype=bool))
                    selected = np.concatenate(selected)[:self.imax]
                    assert(len(selected) == len(newdata[:,0]))
                elif not (len(newdata[:,0]) == len(selected)):
                    # if simulation is still running, subsequent newdata may
                    # be longer
                    self.imax = min(len(selected), len(newdata[:,0]))
                    selected = selected[:self.imax]
                    newdata = newdata[:self.imax,:]
                # select only unique data
                newdata = newdata[selected,:]

            # set time-height data
            setattr( self, field, newdata[:,2:] )
            self._processed.append(field)
            print('  read',field)

        # set time step arrays (should be identical for all fields)
        self.t = np.array(newdata[:,0])
        self.dt = np.array(newdata[:,1])
        

    def _trim_series_if_needed(self,fields_to_check=None):
        """check for inconsistent array lengths and trim if needed"""
        if fields_to_check is None:
            fields_to_check = self._processed
        for field in fields_to_check:
            try:
                getattr(self,field)
            except AttributeError:
                fields_to_check.remove(field)
        field_lengths = [ getattr(self,field).shape[0] for field in fields_to_check ]
        if np.min(field_lengths) < np.max(field_lengths):
            self.imax = np.min(field_lengths)
            # need to prune arrays
            print('Inconsistent averaging field lengths... is simulation still running?')
            print('  truncated field histories from',np.max(field_lengths),'to',self.imax)
            self.t = self.t[:self.imax]
            self.dt = self.dt[:self.imax]
            for field in fields_to_check:
                setattr(self, field, getattr(self,field)[:self.imax,:])

    #==========================================================================
    #
    # CALCULATIONS
    #
    #==========================================================================

    def calculate_TI(self,
                     heights=[],
                     tavg_window=600.0,
                     dt=1.0,
                     SFS=True,
                     verbose=True):
        """Calculate the turbulence intensity (TI) of the resolved
        fluctuations alone or combined fluctuations (including resolved
        and sub-filter scale, SFS).

        INPUTS
            tavg_window     size of window for moving average [s]
            heights         vertical locations at which to calculate TI [m]
            dt              uniform time interval to which to interpolate [s]
            SFS             set to True to include SFS terms

        CALCULATED QUANTITIES
            tavg            uniformly spaced times at which a moving average was calculated
            TIx_hist        variance in x-dir, *_hist.shape == ( len(tavg), len(heights) )
            TIy_hist        variance in y-dir
            TIz_hist        variance in z-dir
            TIdir_hist      variance resolved to flow direction
            TIxyz_hist      variance assuming homogeneous turbulence, calculated from TKE
            TKE_hist        turbulent kinetic energy (TKE)
        """
        try:
            from scipy.ndimage import uniform_filter
        except ImportError:
            print('Moving average calculation uses scipy.ndimage')
            return

        self.get_vars_if_needed('uu_mean','vv_mean','ww_mean','uv_mean','uw_mean','vw_mean')
        if SFS: self.get_vars_if_needed('R11_mean','R22_mean','R33_mean','R12_mean','R13_mean','R23_mean')

        Nout  = len(heights)
        if Nout==0:
            print('Need to specify output heights')
            return

        self._trim_series_if_needed()

        # setup uniform points for interpolation and averaging windows
        Nt = int(np.ceil((self.t[-1]-self.t[0])/dt))
        tuniform = np.arange(1,Nt+1)*dt + self.t[0]
        Navg    = int(tavg_window/dt)
        Navg_2  = int(Navg/2)
        tavg    = tuniform[Navg_2:-Navg_2+1]
        Ntavg   = len(tavg)
        if verbose:
            print('Interpolating to',Nt,'uniformly-spaced data points')
            print('Moving average window:',tavg_window,'s')

        TIx   = np.zeros((Ntavg,Nout))
        TIy   = np.zeros((Ntavg,Nout))
        TIz   = np.zeros((Ntavg,Nout))
        TIdir = np.zeros((Ntavg,Nout))
        TIxyz = np.zeros((Ntavg,Nout))
        TKE   = np.zeros((Ntavg,Nout))

        for ih,z in enumerate(heights):

            # setup interpolation
            idx = np.nonzero(self.hLevelsCell > heights[ih])[0]
            if len(idx)==0: # requested height too high
                k = len(self.hLevelsCell) # extrapolate from two highest points
            elif idx[0]==0: # request height too low
                k = 1 # extrapolate from two lowest points
            else:
                k = idx[0] # interpolate between k-1, k
            frac = (heights[ih] - self.hLevelsCell[k-1]) / (self.hLevelsCell[k] - self.hLevelsCell[k-1])

            # calculate time-averaged velocity profiles
            # 1. interpolate to requested height for all times
            # 2. interpolate from all times to a time history with uniformly-spaced samples
            # 3. apply uniform_filter to perform moving average over the uniformly-spaced samples
            U_mean_interp = self.U_mean[:,k-1] + frac*(self.U_mean[:,k] - self.U_mean[:,k-1]) # length=len(self.t)
            V_mean_interp = self.V_mean[:,k-1] + frac*(self.V_mean[:,k] - self.V_mean[:,k-1])
            W_mean_interp = self.W_mean[:,k-1] + frac*(self.W_mean[:,k] - self.W_mean[:,k-1])
            U_mean_uniform = np.interp( tuniform, self.t, U_mean_interp ) # length=Nt
            V_mean_uniform = np.interp( tuniform, self.t, V_mean_interp )
            W_mean_uniform = np.interp( tuniform, self.t, W_mean_interp )
            UMeanAvg = uniform_filter( U_mean_uniform, Navg )[Navg_2:-Navg_2+1] # length=Ntavg
            VMeanAvg = uniform_filter( V_mean_uniform, Navg )[Navg_2:-Navg_2+1]
            WMeanAvg = uniform_filter( W_mean_uniform, Navg )[Navg_2:-Navg_2+1]

            # calculate time-averaged variances
            uu_mean_interp = self.uu_mean[:,k-1] + frac*(self.uu_mean[:,k] - self.uu_mean[:,k-1]) # length=len(self.t)
            vv_mean_interp = self.vv_mean[:,k-1] + frac*(self.vv_mean[:,k] - self.vv_mean[:,k-1])
            uv_mean_interp = self.uv_mean[:,k-1] + frac*(self.uv_mean[:,k] - self.uv_mean[:,k-1])
            ww_mean_interp = self.ww_mean[:,k-1] + frac*(self.ww_mean[:,k] - self.ww_mean[:,k-1])
            uu_mean_uniform = np.interp( tuniform, self.t, uu_mean_interp ) # length=Nt
            vv_mean_uniform = np.interp( tuniform, self.t, vv_mean_interp )
            uv_mean_uniform = np.interp( tuniform, self.t, uv_mean_interp )
            ww_mean_uniform = np.interp( tuniform, self.t, ww_mean_interp )
            uuMeanAvg = uniform_filter( uu_mean_uniform, Navg )[Navg_2:-Navg_2+1] # length=Ntavg
            vvMeanAvg = uniform_filter( vv_mean_uniform, Navg )[Navg_2:-Navg_2+1]
            uvMeanAvg = uniform_filter( uv_mean_uniform, Navg )[Navg_2:-Navg_2+1]
            wwMeanAvg = uniform_filter( ww_mean_uniform, Navg )[Navg_2:-Navg_2+1]
            if SFS:
                if verbose: print('Adding SFS component')
                R11_mean_interp = self.R11_mean[:,k-1] + frac*(self.R11_mean[:,k] - self.R11_mean[:,k-1]) # length=len(self.t)
                R22_mean_interp = self.R22_mean[:,k-1] + frac*(self.R22_mean[:,k] - self.R22_mean[:,k-1])
                R12_mean_interp = self.R12_mean[:,k-1] + frac*(self.R12_mean[:,k] - self.R12_mean[:,k-1])
                R33_mean_interp = self.R33_mean[:,k-1] + frac*(self.R33_mean[:,k] - self.R33_mean[:,k-1])
                R11_mean_uniform = np.interp( tuniform, self.t, R11_mean_interp ) #length=Nt
                R22_mean_uniform = np.interp( tuniform, self.t, R22_mean_interp )
                R12_mean_uniform = np.interp( tuniform, self.t, R12_mean_interp )
                R33_mean_uniform = np.interp( tuniform, self.t, R33_mean_interp )
                uuMeanAvg += uniform_filter( R11_mean_uniform, Navg )[Navg_2:-Navg_2+1] # length=Ntavg
                vvMeanAvg += uniform_filter( R22_mean_uniform, Navg )[Navg_2:-Navg_2+1]
                uvMeanAvg += uniform_filter( R12_mean_uniform, Navg )[Navg_2:-Navg_2+1]
                wwMeanAvg += uniform_filter( R33_mean_uniform, Navg )[Navg_2:-Navg_2+1]

            Umag = np.sqrt( UMeanAvg**2 + VMeanAvg**2 + WMeanAvg**2 )
            windDir = np.abs( np.arctan2(VMeanAvg,UMeanAvg) )

            # calculate TKE and TI
            TIx[:,ih] = np.sqrt( uuMeanAvg ) / Umag
            TIy[:,ih] = np.sqrt( vvMeanAvg ) / Umag
            TIz[:,ih] = np.sqrt( wwMeanAvg ) / Umag
            TKE[:,ih] = 0.5*( uuMeanAvg + vvMeanAvg + wwMeanAvg )
            TIxyz[:,ih] = np.sqrt( 2./3.*TKE[:,ih] ) / Umag

            TIdir[:,ih] = uuMeanAvg *   np.cos(windDir)**2 \
                        + uvMeanAvg * 2*np.sin(windDir)*np.cos(windDir) \
                        + vvMeanAvg *   np.sin(windDir)**2
            TIdir[:,ih] = np.sqrt(TIdir[:,ih]) / Umag

        # end loop over heights

        # save attributes, shape==(Ntavg,Nout)
        self.TI_heights = heights
        self.tavg       = tavg
        self.TIx        = TIx
        self.TIy        = TIy
        self.TIz        = TIz
        self.TIdir      = TIdir
        self.TIxyz      = TIxyz
        self.TKE        = TKE
        
    def calculate_TI_profile(self,
                             time=9e9,
                             tavg_window=600.0,
                             dt=1.0,
                             SFS=True,
                             verbose=True):
        """Calculate the directional turbulence intensity (TI) profile
        of the resolved fluctuations alone or combined fluctuations
        (including resolved and sub-filter scale, SFS).

        INPUTS
            tavg_window     size of window for moving average [s]
            dt              uniform time interval to which to interpolate [s]
            SFS             set to True to include SFS terms

        CALCULATED QUANTITIES
            TI_profile      profile of variance resolved to flow direction
            TKE_profile     profile of variance resolved to flow direction
        """
        try:
            from scipy.ndimage import uniform_filter
        except ImportError:
            print('Moving average calculation uses scipy.ndimage')
            return

        self.get_vars_if_needed('uu_mean','vv_mean','ww_mean','uv_mean','uw_mean','vw_mean')
        if SFS: self.get_vars_if_needed('R11_mean','R22_mean','R33_mean','R12_mean','R13_mean','R23_mean')

        self._trim_series_if_needed()

        # setup uniform points for interpolation and averaging windows
        Nt = int(np.ceil((self.t[-1]-self.t[0])/dt))
        tuniform = np.arange(1,Nt+1)*dt + self.t[0]
        Navg    = int(tavg_window/dt)
        Navg_2  = int(Navg/2)
        tavg    = tuniform[Navg_2:-Navg_2+1]
        Ntavg   = len(tavg)
        if verbose:
            print('Interpolating to',Nt,'uniformly-spaced data points')
            print('Moving average window:',tavg_window,'s')

        TI_profile = np.zeros((len(self.hLevelsCell)))
        TKE_profile = np.zeros((len(self.hLevelsCell)))

        for ih,z in enumerate(self.hLevelsCell):

            # calculate time-averaged velocity profiles
            # 1. interpolate from all times to a time history with uniformly-spaced samples
            # 2. apply uniform_filter to perform moving average over the uniformly-spaced samples
            U_mean_uniform = np.interp( tuniform, self.t, self.U_mean[:,ih] ) # length=Nt
            V_mean_uniform = np.interp( tuniform, self.t, self.V_mean[:,ih] )
            W_mean_uniform = np.interp( tuniform, self.t, self.W_mean[:,ih] )
            UMeanAvg = uniform_filter( U_mean_uniform, Navg )[Navg_2:-Navg_2+1] # length=Ntavg
            VMeanAvg = uniform_filter( V_mean_uniform, Navg )[Navg_2:-Navg_2+1]
            WMeanAvg = uniform_filter( W_mean_uniform, Navg )[Navg_2:-Navg_2+1]

            # calculate time-averaged variances
            uu_mean_uniform = np.interp( tuniform, self.t, self.uu_mean[:,ih] ) # length=Nt
            vv_mean_uniform = np.interp( tuniform, self.t, self.vv_mean[:,ih] )
            uv_mean_uniform = np.interp( tuniform, self.t, self.uv_mean[:,ih] )
            ww_mean_uniform = np.interp( tuniform, self.t, self.ww_mean[:,ih] )
            uuMeanAvg = uniform_filter( uu_mean_uniform, Navg )[Navg_2:-Navg_2+1] # length=Ntavg
            vvMeanAvg = uniform_filter( vv_mean_uniform, Navg )[Navg_2:-Navg_2+1]
            uvMeanAvg = uniform_filter( uv_mean_uniform, Navg )[Navg_2:-Navg_2+1]
            wwMeanAvg = uniform_filter( ww_mean_uniform, Navg )[Navg_2:-Navg_2+1]
            if SFS:
                R11_mean_uniform = np.interp( tuniform, self.t, self.R11_mean[:,ih] ) #length=Nt
                R22_mean_uniform = np.interp( tuniform, self.t, self.R22_mean[:,ih] )
                R12_mean_uniform = np.interp( tuniform, self.t, self.R12_mean[:,ih] )
                R33_mean_uniform = np.interp( tuniform, self.t, self.R33_mean[:,ih] )
                uuMeanAvg += uniform_filter( R11_mean_uniform, Navg )[Navg_2:-Navg_2+1] # length=Ntavg
                vvMeanAvg += uniform_filter( R22_mean_uniform, Navg )[Navg_2:-Navg_2+1]
                uvMeanAvg += uniform_filter( R12_mean_uniform, Navg )[Navg_2:-Navg_2+1]
                wwMeanAvg += uniform_filter( R33_mean_uniform, Navg )[Navg_2:-Navg_2+1]

            Umag = np.sqrt( UMeanAvg**2 + VMeanAvg**2 + WMeanAvg**2 )
            windDir = np.abs( np.arctan2(VMeanAvg,UMeanAvg) )

            # calculate TKE and TI
            #TIx_profile[ih] = np.sqrt( uuMeanAvg ) / Umag
            #TIy_profile[ih] = np.sqrt( vvMeanAvg ) / Umag
            #TIz_profile[ih] = np.sqrt( wwMeanAvg ) / Umag
            #TIxyz_profile[ih] = np.sqrt( 2./3.*TKE[:,ih] ) / Umag

            itime = np.argmin(np.abs(time - self.tavg))
            TKE_profile[ih] = 0.5*( uuMeanAvg[itime] + vvMeanAvg[itime] + wwMeanAvg[itime] )
            TI_profile[ih] = uuMeanAvg[itime] *   np.cos(windDir[itime])**2 \
                           + uvMeanAvg[itime] * 2*np.sin(windDir[itime])*np.cos(windDir[itime]) \
                           + vvMeanAvg[itime] *   np.sin(windDir[itime])**2
            TI_profile[ih] = np.sqrt(TI_profile[ih]) / Umag[itime]

        # end loop over heights

        # save attributes, shape==(Ntavg,Nout)
        self.tavg_profile = tavg
        self.TI_profile   = TI_profile
        self.TKE_profile  = TKE_profile
        
    def calculate_shear(self,
                        heights=[20.0,40.0,80.0],
                        zref=80.0,
                        Uref=8.0,
                        interp=None,
                        verbose=True):
        """Estimate the shear from the average streamwise velocity from
        the final time step. Sets self.approxWindProfile to the fitted
        wind profile.

        INPUTS
            heights     list of points to use to fit the wind profile
            zref        power-law reference height
            Uref        power-law reference velocity
            interp      None, or 'kind' input to scipy.interpolate.interp1d

        OUTPUTS
            alpha       power law wind profile exponent
        """
        from scipy.interpolate import interp1d
        Uh = np.sqrt( self.U_mean**2 + self.V_mean**2 )[-1,:] # horizontal wind

        if interp is not None:
            Ufn = interp1d( self.hLevelsCell, Uh, kind=interp )
            U = Ufn(heights)
            self.approxHeights = heights
            self.approxU = U # at heights
            self.approxUfn = Ufn
        else:
            # find nearest cell without interpolating
            idx = [ np.argmin(np.abs(h-self.hLevelsCell)) for h in sorted(heights) ]
            assert(len(set(idx)) >= 2)
            heights = self.hLevelsCell[idx]
            U = Uh[idx]

        if verbose:
            print('Estimating shear coefficient for Uref=',Uref,'and zref=',zref,':')
            print('     U=',U,'m/s')
            print('  at z=',heights,'m')
        lnz = np.log( np.array(heights)/zref )
        lnU = np.log( U/Uref )
        alpha = lnz.dot(lnU) / lnz.dot(lnz)

        self.Uh = Uh
        self.approxWindProfile = Uref * (self.hLevelsCell/zref)**alpha
        self.alpha = alpha

        return self.alpha

    def calculate_veer(self, zhub=80.0, D=126.0, verbose=True):
        """Estimate the veer from the average streamwise velocity from
        the final time step. Also calculates the height-varying wind
        direction, saved in self.windDir [deg].

        INPUTS
            zhub        hub height [m]
            D           rotor diameter [m]

        OUTPUTS
            veer        veer angle [deg]
                        >0 implies clockwise change in wind direction seen from above
        """
        dir = np.arctan2( -self.U_mean, -self.V_mean )[-1,:]
        dir[dir<0] += 2*np.pi
        dir *= 180.0/np.pi

        self.windDir = dir
        
        idxs = (self.hLevelsCell >= zhub-D/2) & (self.hLevelsCell <= zhub+D/2)
        if verbose:
            z = self.hLevelsCell[idxs]
            #for zi,angi,tf in zip(self.hLevelsCell,dir,idxs):
            #    print(zi,'m ',angi,'deg ',tf)
            print('Estimating shear between z=',z[0],'and',z[-1],'m')
        rotorWindDir = dir[idxs]
        self.veer = rotorWindDir[-1] - rotorWindDir[0]

        return self.veer

    def calculate_Tgrad(self, zi):
        """Calculate the temperature gradient at the inversion height
        and at the top of the computational domain.

        INPUTS
            zi          inversion height [m]

        OUTPUTS
            dTdz_inv    temperature gradient at the inversion [deg K/m]
            dTdz_upper  upper temperature gradient [deg K/m]
        """
        self.get_vars_if_needed('T_mean')

        ii = np.argmin(np.abs(self.hLevelsCell-zi))
        dTdz_inv = (self.T_mean[-1,ii+1] - self.T_mean[-1,ii]) \
                 / (self.hLevelsCell[ii+1] - self.hLevelsCell[ii])

        dTdz_upper = (self.T_mean[-1,-1] - self.T_mean[-1,-2]) \
                   / (self.hLevelsCell[-1] - self.hLevelsCell[-2])

        return dTdz_inv, dTdz_upper

    def calculate_Richardson(self, g=9.81, zref=90.0, D=126.0, verbose=True):
        """Estimate the Richardson number from the averaged T profile
        at the final time step which should range between -O(0.01) to
        +O(0.1), for unstable to stable. Difference formula used are
        second-order accurate.

        INPUTS
            g           gravity [m/s^2]
            zref        reference height used to locate the top of the rotor [m]
            D           rotor diameter used to locate the top of the rotor [m]

        OUTPUTS
            Ri          Richardson number
        """
        self.get_vars_if_needed('T_mean')

        z = self.hLevelsCell
        T = self.T_mean[-1,:]
        U = np.sqrt( self.U_mean[-1,:]**2 + self.V_mean[-1,:]**2 )
        spacings = z[1:3] - z[0:2]
        assert( spacings[0] == spacings[1] )
        dz = spacings[0]

        rotorTop = zref + D/2
        idxs = z<=rotorTop
        Tmean = np.mean( T[idxs] )

        centralFormula = np.array([-1,0,1])/(2*dz)
        dTdz1 = centralFormula.dot(T[:3])
        dUdz1 = centralFormula.dot(U[:3])

        oneSidedFormula = np.array([-3,4,-1])/(2*dz)
        dTdz0 = oneSidedFormula.dot(T[:3])
        dUdz0 = oneSidedFormula.dot(U[:3])

        Tsurf = (T[1]-T[0])/(z[1]-z[0])*(-z[0]) + T[0]

        # extrapolate derivatives to the ground
        dTdz = (dTdz1-dTdz0)/dz * (-z[0]) + dTdz0
        dUdz = (dUdz1-dUdz0)/dz * (-z[0]) + dUdz0

        if verbose:
            print('Calculating Ri with:')
            print('  mean T :',Tmean)
            print('  dT/dz at z=',z[1],':',dTdz1,' (finite difference)')
            print('  dT/dz at z=',z[0],':',dTdz0,' (finite difference)')
            print('  dT/dz at z=0:',dTdz,' (extrapolated)')
            print('  dU/dz at z=',z[1],':',dUdz1,' (finite difference)')
            print('  dU/dz at z=',z[0],':',dUdz0,' (finite difference)')
            print('  dU/dz at z=0:',dUdz,' (extrapolated)')

        self.Tsurf = Tsurf
        self.dUdz_surf = dUdz
        self.dTdz_surf = dTdz
        self.Ri = g/Tmean * dTdz / dUdz**2

        return self.Ri

    def rotate_tensors(self,itime=None):
        """Rotate the resolved and SFS stress tensors, using tensor
        transformation rules, at the specified time index (otherwise
        all times are calculated--this is slow).
        """
        if not self.rotated:
            self.uu_meanR = np.zeros(self.uu_mean.shape)
            self.vv_meanR = np.zeros(self.vv_mean.shape)
            self.ww_meanR = np.zeros(self.ww_mean.shape)
            self.uv_meanR = np.zeros(self.uv_mean.shape)
            self.uw_meanR = np.zeros(self.uw_mean.shape)
            self.vw_meanR = np.zeros(self.vw_mean.shape)
            self.R11_meanR = np.zeros(self.R11_mean.shape)
            self.R22_meanR = np.zeros(self.R22_mean.shape)
            self.R33_meanR = np.zeros(self.R33_mean.shape)
            self.R12_meanR = np.zeros(self.R12_mean.shape)
            self.R13_meanR = np.zeros(self.R13_mean.shape)
            self.R23_meanR = np.zeros(self.R23_mean.shape)

        if itime is None:
            print('Rotating stress tensors for all times...')
            tindices = range(len(self.t))
        else:
            print('Rotating stress tensors at t={:f} s'.format(self.t[itime]))
            tindices = [itime]

        ang = np.arctan2(self.V_mean, self.U_mean) # NOT compass wind dir

        for iz in range(len(self.hLevelsCell)):
            for it in tindices:
                # resolved stress tensor
                S = np.array([[ self.uu_mean[it,iz],  self.uv_mean[it,iz],  self.uw_mean[it,iz]],
                              [ self.uv_mean[it,iz],  self.vv_mean[it,iz],  self.vw_mean[it,iz]],
                              [ self.uw_mean[it,iz],  self.vw_mean[it,iz],  self.ww_mean[it,iz]]])
                # Reynolds stress tensor
                R = np.array([[self.R11_mean[it,iz], self.R12_mean[it,iz], self.R13_mean[it,iz]],
                              [self.R12_mean[it,iz], self.R22_mean[it,iz], self.R23_mean[it,iz]],
                              [self.R13_mean[it,iz], self.R23_mean[it,iz], self.R33_mean[it,iz]]])

                S = rotate_tensor_z(S, ang[it,iz])
                R = rotate_tensor_z(R, ang[it,iz])

                self.uu_meanR[it,iz] = S[0,0]
                self.uv_meanR[it,iz] = S[0,1]
                self.uw_meanR[it,iz] = S[0,2]
                self.vv_meanR[it,iz] = S[1,1]
                self.vw_meanR[it,iz] = S[1,2]
                self.ww_meanR[it,iz] = S[2,2]
                self.R11_meanR[it,iz] = R[0,0]
                self.R12_meanR[it,iz] = R[0,1]
                self.R13_meanR[it,iz] = R[0,2]
                self.R22_meanR[it,iz] = R[1,1]
                self.R23_meanR[it,iz] = R[1,2]
                self.R33_meanR[it,iz] = R[2,2]

        self.rotated = True

    #==========================================================================
    #
    # PLOTS
    #
    #==========================================================================

    def plot_TI_history(self,ax=None,savefig=None):
        """Plots TI history at all height at which TI was calculated
        An optional image name may be specified as 'savefig'
        """
        if ax is None:
            fig,ax = plt.subplots()
        else:
            fig = plt.gcf()
        cmap = cm.get_cmap('viridis')
        for ih,z in enumerate(self.TI_heights):
            f = (z-self.TI_heights[0]) / self.TI_heights[-1]
            ax.plot(self.tavg, 100.0*self.TIdir[:,ih],
                    color=cmap(f),
                    label='z={:.1f} m'.format(z))
        ax.set_xlabel(r'Time [s]')
        ax.set_ylabel(r'Turbulence Intensity [%]')
        ax.legend(loc='best',fontsize='small')
        if savefig is not None:
            fig.savefig(savefig,bbox_inches='tight')
        return fig, ax

    def plot_TKE_history(self,ax=None,savefig=None):
        """Plots TKE history at all height at which TI was calculated
        An optional image name may be specified as 'savefig'
        """
        if ax is None:
            fig,ax = plt.subplots()
        else:
            fig = plt.gcf()
        cmap = cm.get_cmap('viridis')
        for ih,z in enumerate(self.TI_heights):
            f = (z-self.TI_heights[0]) / self.TI_heights[-1]
            ax.plot(self.tavg, self.TKE[:,ih],
                    color=cmap(f),
                    label='z={:.1f} m'.format(z))
        ax.set_xlabel(r'Time [s]')
        ax.set_ylabel(r'Turbulent Kinetic Energy [m$^2$/s$^2$]')
        ax.legend(loc='best',fontsize='small')
        if savefig is not None:
            fig.savefig(savefig,bbox_inches='tight')
        return fig, ax

    def plot_UVW_profile(self,time=9e9,ax=None,savefig=None,**kwargs):
        """Plot profiles of wind velocity components at the time instant
        nearest to the specified time
        """
        if ax is None:
            fig,ax = plt.subplots()
        else:
            fig = plt.gcf()
        itime = np.argmin(np.abs(time - self.t))
        ax.plot(self.U_mean[itime,:], self.hLevelsCell, label=r'$U$', **kwargs)
        ax.plot(self.V_mean[itime,:], self.hLevelsCell, label=r'$V$', **kwargs)
        ax.plot(self.W_mean[itime,:], self.hLevelsCell, label=r'$W$', **kwargs)
        ax.set_xlabel(r'Velocity [m/s]')
        ax.set_ylabel(r'Height [m]')
        ax.legend(loc='best')
        if savefig is not None:
            fig.savefig(savefig,bbox_inches='tight')
        return fig, ax

    def plot_windspeed_profile(self,time=9e9,ax=None,savefig=None,**kwargs):
        """Plot profiles of wind speed magnitude at the time instant
        nearest to the specified time
        """
        if ax is None:
            fig,ax = plt.subplots()
        else:
            fig = plt.gcf()
        itime = np.argmin(np.abs(time - self.t))
        windMag = np.sqrt( self.U_mean[itime,:]**2 + self.V_mean[itime,:]**2 )
        ax.plot(windMag, self.hLevelsCell, **kwargs)
        ax.set_xlabel(r'Horizontal Velocity [m/s]')
        ax.set_ylabel(r'Height [m]')
        if savefig is not None:
            fig.savefig(savefig,bbox_inches='tight')
        return fig, ax

    def plot_winddirection_profile(self,time=9e9,ax=None,savefig=None,**kwargs):
        """Plot profiles of wind direction at the time instant nearest
        to the specified time
        """
        if ax is None:
            fig,ax = plt.subplots()
        else:
            fig = plt.gcf()
        itime = np.argmin(np.abs(time - self.t))
        windDir = 180.0/np.pi * np.arctan2(-self.U_mean[itime,:], -self.V_mean[itime,:])
        meanWindDir = np.mean(windDir)
        if meanWindDir < 0:
            meanWindDir += 360.0
            windDir += 360.0
        ax.plot(windDir, self.hLevelsCell, **kwargs)
        cur_xlim = ax.get_xlim()
        xlim = [ min(cur_xlim[0],np.round(meanWindDir-1.0)),
                 max(cur_xlim[1],np.round(meanWindDir+1.0)) ]
        ax.set_xlim(xlim)
        ax.set_xlabel(r'Wind Direction [deg]')
        ax.set_ylabel(r'Height [m]')
        if savefig is not None:
            fig.savefig(savefig,bbox_inches='tight')
        return fig, ax

    def plot_T_profile(self,time=9e9,ax=None,savefig=None,**kwargs):
        """Plot temperature profile at the time instant nearest to the
        specified time
        """
        if ax is None:
            fig,ax = plt.subplots()
        else:
            fig = plt.gcf()
        itime = np.argmin(np.abs(time - self.t))
        ax.plot(self.T_mean[itime,:], self.hLevelsCell, label=r'$T$', **kwargs)
        ax.set_xlabel(r'Temperature [K]')
        ax.set_ylabel(r'Height [m]')
        if savefig is not None:
            fig.savefig(savefig,bbox_inches='tight')
        return fig, ax

    def plot_variance_profile(self,time=9e9,rotated=False,ax=None,savefig=None,**kwargs):
        """Plot profiles of variance at the time instant nearest to the
        specified time
        """
        if ax is None:
            fig,ax = plt.subplots()
        else:
            fig = plt.gcf()
        itime = np.argmin(np.abs(time - self.t))
        self.get_vars_if_needed('uu_mean','vv_mean','ww_mean')
        if rotated:
            self.rotate_tensors(itime)
            ax.plot(self.uu_meanR[itime,:], self.hLevelsCell, label=r"$<u'u'>$", **kwargs)
            ax.plot(self.vv_meanR[itime,:], self.hLevelsCell, label=r"$<v'v'>$", **kwargs)
            ax.plot(self.ww_meanR[itime,:], self.hLevelsCell, label=r"$<w'w'>$", **kwargs)
            ax.set_xlabel(r'Variance (wind-aligned) [m$^2$/s$^2$]')
        else:
            ax.plot(self.uu_mean[itime,:], self.hLevelsCell, label=r"$<u'u'>$", **kwargs)
            ax.plot(self.vv_mean[itime,:], self.hLevelsCell, label=r"$<v'v'>$", **kwargs)
            ax.plot(self.ww_mean[itime,:], self.hLevelsCell, label=r"$<w'w'>$", **kwargs)
            ax.set_xlabel(r'Variance [m$^2$/s$^2$]')
        ax.set_ylabel(r'Height [m]')
        ax.legend(loc='best')
        if savefig is not None:
            fig.savefig(savefig,bbox_inches='tight')
        return fig, ax

    def plot_covariance_profile(self,time=9e9,rotated=False,ax=None,savefig=None,**kwargs):
        """Plot profiles of covariance at the time instant nearest to the
        specified time
        """
        if ax is None:
            fig,ax = plt.subplots()
        else:
            fig = plt.gcf()
        itime = np.argmin(np.abs(time - self.t))
        self.get_vars_if_needed('uv_mean','uw_mean','vw_mean','Tw_mean')
        if rotated:
            self.rotate_tensors(itime)
            ax.plot(self.uv_meanR[itime,:], self.hLevelsCell, label=r"$<u'v'>$", **kwargs)
            ax.plot(self.uw_meanR[itime,:], self.hLevelsCell, label=r"$<u'w'>$", **kwargs)
            ax.plot(self.vw_meanR[itime,:], self.hLevelsCell, label=r"$<v'w'>$", **kwargs)
            ax.plot(self.Tw_mean[itime,:], self.hLevelsCell, label=r"$<T'w'>$", **kwargs)
            ax.set_xlabel(r'Covariance (wind-aligned) [m$^2$/s$^2$], [K-m/s]')
        else:
            ax.plot(self.uv_mean[itime,:], self.hLevelsCell, label=r"$<u'v'>$", **kwargs)
            ax.plot(self.uw_mean[itime,:], self.hLevelsCell, label=r"$<u'w'>$", **kwargs)
            ax.plot(self.vw_mean[itime,:], self.hLevelsCell, label=r"$<v'w'>$", **kwargs)
            ax.plot(self.Tw_mean[itime,:], self.hLevelsCell, label=r"$<T'w'>$", **kwargs)
            ax.set_xlabel(r'Covariance [m$^2$/s$^2$], [K-m/s]')
        ax.set_ylabel(r'Height [m]')
        ax.legend(loc='best')
        if savefig is not None:
            fig.savefig(savefig,bbox_inches='tight')
        return fig, ax

    def plot_SFS_normalstress_profile(self,time=9e9,rotated=False,ax=None,savefig=None,**kwargs):
        """Plot profiles of the modeled sub-filter normal stresses at
        the time instant nearest to the specified time
        """
        if ax is None:
            fig,ax = plt.subplots()
        else:
            fig = plt.gcf()
        itime = np.argmin(np.abs(time - self.t))
        self.get_vars_if_needed('R11_mean','R22_mean','R33_mean')
        if rotated:
            self.rotate_tensors(itime)
            ax.plot(self.R11_meanR[itime,:], self.hLevelsCell, label=r"$R_{11}$", **kwargs)
            ax.plot(self.R22_meanR[itime,:], self.hLevelsCell, label=r"$R_{22}$", **kwargs)
            ax.plot(self.R33_meanR[itime,:], self.hLevelsCell, label=r"$R_{33}$", **kwargs)
            ax.set_xlabel(r'SFS Normal Stresses (wind-aligned) [m$^2$/s$^2$]')
        else:
            ax.plot(self.R11_mean[itime,:], self.hLevelsCell, label=r"$R_{11}$", **kwargs)
            ax.plot(self.R22_mean[itime,:], self.hLevelsCell, label=r"$R_{22}$", **kwargs)
            ax.plot(self.R33_mean[itime,:], self.hLevelsCell, label=r"$R_{33}$", **kwargs)
            ax.set_xlabel(r'SFS Normal Stresses [m$^2$/s$^2$]')
        ax.set_ylabel(r'Height [m]')
        ax.legend(loc='best')
        if savefig is not None:
            fig.savefig(savefig,bbox_inches='tight')
        return fig, ax

    def plot_SFS_shearstress_profile(self,time=9e9,rotated=False,ax=None,savefig=None,**kwargs):
        """Plot profiles of the modeled sub-filter shear stresses at the
        time instant nearest to the specified time
        """
        if ax is None:
            fig,ax = plt.subplots()
        else:
            fig = plt.gcf()
        itime = np.argmin(np.abs(time - self.t))
        self.get_vars_if_needed('R12_mean','R13_mean','R23_mean')
        if rotated:
            self.rotate_tensors(itime)
            ax.plot(self.R12_meanR[itime,:], self.hLevelsCell, label=r"$R_{12}$", **kwargs)
            ax.plot(self.R13_meanR[itime,:], self.hLevelsCell, label=r"$R_{13}$", **kwargs)
            ax.plot(self.R23_meanR[itime,:], self.hLevelsCell, label=r"$R_{23}$", **kwargs)
            ax.set_xlabel(r'SFS Shear Stresses (wind-aligned) [m$^2$/s$^2$]')
        else:
            ax.plot(self.R12_mean[itime,:], self.hLevelsCell, label=r"$R_{12}$", **kwargs)
            ax.plot(self.R13_mean[itime,:], self.hLevelsCell, label=r"$R_{13}$", **kwargs)
            ax.plot(self.R23_mean[itime,:], self.hLevelsCell, label=r"$R_{23}$", **kwargs)
            ax.set_xlabel(r'SFS Shear Stresses [m$^2$/s$^2$]')
        ax.set_ylabel(r'Height [m]')
        ax.legend(loc='best')
        if savefig is not None:
            fig.savefig(savefig,bbox_inches='tight')
        return fig, ax

    #==========================================================================
    #
    # DATA I/O
    #
    #==========================================================================

    def save_TI_history(self,prefix='TIhist',writeTKE=True):
        """Writes out one csv file per height at which TI was calculated"""
        for ih,z in enumerate(self.TI_heights):
            fname = '{:s}_z{:.1f}.csv'.format(prefix,z)
            if writeTKE:
                header = 'Time,TI,TKE'
                TIdata = (self.tavg, self.TIdir[:,ih], self.TKE[:,ih])
            else:
                header = 'Time,TI'
                TIdata = (self.tavg, self.TIdir[:,ih])
            try:
                np.savetxt(fname,
                           np.vstack(TIdata).T,
                           delimiter=',',
                           header=header)
                print('wrote',fname)
            except IOError:
                print('Error:',fname,'could not be written')

    def save_profile(self,time=9e9,fname='averagingProfiles.csv'):
        """Writes out a csv file with the planar-averaged profile
        nearest to the specified instance in time. 
        """
        itime = np.argmin(np.abs(time - self.t))
        print('Outputting averaged profile at',self.t[itime])
        try:
            np.savetxt( fname,
                    np.vstack((
                        self.hLevelsCell,
                        self.U_mean[itime,:],
                        self.V_mean[itime,:],
                        self.W_mean[itime,:],
                        self.T_mean[itime,:],
                        self.uu_mean[itime,:],
                        self.vv_mean[itime,:],
                        self.ww_mean[itime,:],
                        self.uv_mean[itime,:],
                        self.uw_mean[itime,:],
                        self.vw_mean[itime,:],
                        self.Tw_mean[itime,:],
                        self.R11_mean[itime,:],
                        self.R22_mean[itime,:],
                        self.R33_mean[itime,:],
                        self.R12_mean[itime,:],
                        self.R13_mean[itime,:],
                        self.R23_mean[itime,:],
                        )).T,
                    header='z,U,V,W,T,uu,vv,ww,uv,uw,vw,Tw,R11,R22,R33,R12,R13,R23',
                    delimiter=',' )
            print('wrote',fname)
        except IOError:
            print('Error:',fname,'could not be written')

    def to_csv(self,fname,**kwargs):
        """Write out specified range of times in a pandas dataframe

        kwargs: see PlanarAverages.to_pandas()
        """
        df = self.to_pandas(**kwargs)
        print('Dumping dataframe to',fname)
        df.to_csv(fname)

    def to_pandas(self,itime=None,fields=None,dtype=None):
        """Create pandas dataframe for the specified range of times

        Inputs
        ------
        itime: integer, list
            Time indice(s) to write out; if None, all times are output
        fields: list, dict
            Name of field variables to write out; if None, all variables
            that have been processed are written out (including
            variables that have been implicitly read through calls to 
            calculate_* routines). If a dictionary is provided, the keys
            will be the output column names.
        dtype: type
            Single datatype to which to cast all fields
        """
        import pandas as pd
        # output all vars
        if (fields is not None) and (fields.lower() == 'all'):
            print('All fields requested')
            self.get_vars_if_needed(*all_vars)
            self._trim_series_if_needed()
        # select time range
        if itime is None:
            tindices = range(len(self.t))
        else:
            try:
                iter(itime)
            except TypeError:
                # specified single time index
                tindices = [itime]
            else:
                # specified list of indices
                tindices = itime
        # create dataframes for each time (with height as secondary index)
        print('Creating dataframe for',self.t[tindices])
        dflist = []
        for i in tindices:
            if hasattr(fields, '__iter__') and  not isinstance(fields, str):
                # have an iterable object, write out specified fields
                data = { var: getattr(self,var)[i,:] for var in fields }
            elif isinstance(fields, dict):
                # write out specified fields with custom column names
                data = { col: getattr(self,var)[i,:] for col,var in fields.items() }
            else:
                # write out all fields that have been processed
                data = { var: getattr(self,var)[i,:] for var in self._processed }
            data['z'] = self.hLevelsCell
            df = pd.DataFrame(data=data,dtype=dtype)
            df['t'] = self.t[i]
            dflist.append(df)
        df = pd.concat(dflist)
        df = df.set_index(['t','z'])
        return df

"""end of class PlanarAverages"""


def rotate_tensor_z(S,angle):
    """Rotate tensor S about z axis by angle [rad]"""
    R = np.array([[ np.cos(angle), np.sin(angle), 0.],
                  [-np.sin(angle), np.cos(angle), 0.],
                  [            0.,            0., 1.]])
    return np.matmul(np.matmul(R,S), R.T)

