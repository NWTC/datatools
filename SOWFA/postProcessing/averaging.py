#!/usr/bin/python
"""
For processing SOWFA cell averages in postProcessing/averaging
Based on original SOWFA-Tools in MATLAB

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
import os
import numpy as np

def read(*args,**kwargs):
    """ Specify variables to read with keyword 'varList'
    """
    return averagingData(*args,**kwargs)

allAveragingVars = [
        'U_mean','V_mean','W_mean','T_mean',
        'uu_mean', 'vv_mean', 'ww_mean', 'uv_mean', 'uw_mean', 'vw_mean',
        'R11_mean','R22_mean','R33_mean','R12_mean','R13_mean','R23_mean',
        ]

class averagingData(object):

    def __init__(self,*args,**kwargs):
        """ Find and process all time directories
        """# {{{
        self._processed = []
        self.hLevelsCell = None
        self.simTimeDirs = [] # output time names
        self.simStartTimes = [] # start or restart simulation times
        
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
                print 'Checking directory',arg#,'with',listing
                # specified a directory containing output (time) subdirectories
                for dirname in listing:
                    if not os.path.isdir(arg+os.sep+dirname): continue
                    #print '  checking subdirectory',dirname
                    try:
                        startTime = float(dirname)
                        if 'hLevelsCell' in os.listdir(arg+os.sep+dirname):
                            self.simTimeDirs.append( arg+os.sep+dirname )
                            self.simStartTimes.append( startTime )
                    except ValueError: pass # dirname is not a number

        # sort results
        self.simTimeDirs = [ x[1] for x in sorted(zip(self.simStartTimes,self.simTimeDirs)) ]
        self.simStartTimes.sort()

        print 'Simulation (re)start times:',self.simStartTimes

        # process all output dirs
        #for idir,tdir in enumerate( self.simTimeDirs ):
        #    print 'Processing',tdir
        #    self._process(tdir)
        if len(self.simTimeDirs) > 0:
            self._processDirs( self.simTimeDirs, **kwargs )
        else:
            print 'No averaging time directories found!'
    # }}}

    def __repr__(self):
        s = 'SOWFA postProcessing: averaging data'
        for t,d in zip(self.simStartTimes,self.simTimeDirs):
            fullpath = os.path.realpath(os.curdir) + os.sep + d
            s += '\n  {:f}\t{:s}'.format(t,fullpath)
        return s

    def getVarsIfNeeded(self,*args,**kwargs): # reread=False (default keywords after *args not possible in Python 2)
        """ Read in specified list of variables
        """
        varList = []
        if 'reread' in kwargs and kwargs['reread']:
            for var in args:
                if var in self._processed: self._processed.remove(var)
                varList.append(var)
            print 'Rereading variables',varList
        else:
            for var in args:
                if var not in self._processed: varList.append(var)
        if len(varList)==0: return

        self._processDirs( self.simTimeDirs, varList=varList )

    ### DEPRECATED ###
    # This reads each field from one time directory (specified) at a time; several times slower than _processDirs
    def _process(self,tdir):
        """ Reads all files within an averaging output time directory,
        presumably containing hLevelsCell and other cell-averaged
        quantities. An object attribute corresponding to the averaged
        output name is updated, e.g.:
            ${timeDir}/U_mean is appended to the array self.U_mean

        Typically, objects have shape (Nt,Nz).
        """# {{{
        allOutputs = os.listdir(tdir)
        outputs = []
        for field in allOutputs:
            if field=='hLevelsCell': continue
            if field in allAveragingVars:
                outputs.append( out )

        # process hLevelsCell first, verify we have the same cells
        with open(tdir+os.sep+'hLevelsCell','r') as f:
            line = f.readline()
        z = np.array([ float(val) for val in line.split() ])
        if getattr(self,'hLevelsCell') is None: # this is the first time dir
            self.hLevelsCell = z
        elif not np.all( self.hLevelsCell == z ):
            print 'Error: cell levels do not match'
            return

        # check that we have the same amount of data
        Nlines = []
        for field in outputs:
            output = tdir + os.sep + field
            if not os.path.isfile(output): continue

            with open(output,'r') as f:
                for i,line in enumerate(f): pass
                Nlines.append(i+1)
                line = line.split() # check final line for the right number of values
                if not len(line) == len(z)+2: # t,dt,f_1,f_2,...,f_N for N heights
                    print 'z',z
                    print 'line',line
                    print 'Error: number of output points inconsistent with hLevelsCell in',output
                    return

        if not np.min(Nlines) == np.max(Nlines):
            print 'Error: number of output times do not match in all files'
            return
        N = Nlines[0]

        # NOW process all data
        for field in outputs:
            if field in self._processed: continue
            else: self._processed.append(field)

            output = tdir + os.sep + field
            if not os.path.isfile(output): continue

            with open(output,'r') as f:
                lines = f.readlines()
            newdata = np.array([ [ float(val) for val in line.split()] for line in lines ]) # newdata.shape = (Nt,Nz+2)

            try:
                olddata = getattr(self,field)
                setattr( self, field, np.concatenate((olddata,newdata[:,2:])) )
                print '  concatenated',field
            except AttributeError:
                setattr( self, field, newdata[:,2:] )
                print '  created',field

        # add additional times
        try:
            self.t = np.concatenate((self.t, newdata[:,0]))
            self.dt = np.concatenate((self.dt, newdata[:,1]))
        except AttributeError: # this is the first time dir
            self.t = np.array( newdata[:,0] )
            self.dt = np.array( newdata[:,1] )
    # }}}

    # This reads each field from all time directories simultaneously; several times faster than _process!
    def _processDirs(self,tdirList,varList=['U_mean','V_mean','W_mean','T_mean']):
        """ Reads all files within an averaging output time directory,
        presumably containing hLevelsCell and other cell-averaged
        quantities. An object attribute corresponding to the averaged
        output name is updated, e.g.:
            ${timeDir}/U_mean is appended to the array self.U_mean

        Typically, objects have shape (Nt,Nz).
        """# {{{
        outputs = []
        if isinstance( varList, (str,unicode) ):
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

        # check that we have the same amount of data
        for tdir in tdirList:
            Nlines = []
            for field in outputs:
                output = tdir + os.sep + field
                if not os.path.isfile(output):
                    print 'Error:',output,'not found'
                    return

                with open(output,'r') as f:
                    for i,line in enumerate(f): pass
                    Nlines.append(i+1)
                    line = line.split() # check final line for the right number of values
                    if not len(line) == len(self.hLevelsCell)+2: # t,dt,f_1,f_2,...,f_N for N heights
                        print 'z',z
                        print 'line',line
                        print 'Error: number of output points inconsistent with hLevelsCell in',output
                        return
            if not np.min(Nlines) == np.max(Nlines):
                print 'Warning: number of output times do not match in all files'
        N = Nlines[0]

        # NOW process all data
        for field in outputs:
            arrays = [ np.loadtxt( tdir + os.sep + field ) for tdir in tdirList ]
            newdata = np.concatenate(arrays)
            setattr( self, field, newdata[:,2:] )

            self._processed.append(field)
            print '  read',field

        self.t = np.array( newdata[:,0] )
        self.dt = np.array( newdata[:,1] )
        # }}}
        return None

    def calcTI(self,
            heights=[],
            avg_time=None,
            avg_width=300,
            SFS=True,
            verbose=True):
        """ Calculate the turbulence intensity (TI) of the resolved
        fluctuations alone or combined fluctuations (including resolved
        and sub-filter scale, SFS). 
        
        Based on ABLTools/variances_avg_cell.m

        INPUTS
            avg_time        time at which to calculate TI [s]
            avg_width       half of the number of samples over which to perform the statistics
            heights         vertical locations at which to calculate TI [m]
            SFS             set to True to include SFS terms

        OUTPUTS
            TIx             variance in x-dir
            TIy             variance in y-dir
            TIz             variance in z-dir
            TIdir           variance resolved to flow direction
            TIxyz           variance assuming homogeneous turbulence, calculated from TKE
            TKE             turbulent kinetic energy (TKE)
        """# {{{
        self.getVarsIfNeeded('uu_mean','vv_mean','ww_mean','uv_mean','uw_mean','vw_mean')
        if SFS: self.getVarsIfNeeded('R11_mean','R22_mean','R33_mean','R12_mean','R13_mean','R23_mean')

        try:
            Nout = len(heights)
        except TypeError: # specified float instead of list of floats
            heights = [heights]
            Nout = 1
        if Nout==0:
            print 'Need to specify output heights'
            return
        TIx   = np.zeros(Nout)
        TIy   = np.zeros(Nout)
        TIz   = np.zeros(Nout)
        TIdir = np.zeros(Nout)
        TIxyz = np.zeros(Nout)
        TKE   = np.zeros(Nout)

        if avg_time is None:
            avg_time = self.t[-1]
        if verbose: print 'Average at time',avg_time

        dtmin = np.min(self.dt)
        dtmax = np.max(self.dt)
        if dtmin==dtmax:
            if verbose: print 'Constant dt, averaging window :',2*avg_width*dtmin,'s'
        else:
            dtmean = np.mean(self.dt)
            if verbose: 
                print 'Variable dt, approximate averaging window :',avg_width*dtmean,'s'
                print '  dt min/mean/max=',dtmin,dtmean,dtmax

        i = np.argmin(np.abs(self.t-avg_time))
        js = max( i-avg_width  , 0 )
        je = min( i+avg_width+1, len(self.t) )
        dtsub = self.dt[js:je]
        if verbose: print 'Processing range js,je,sum(dt):',js,je,np.sum(dtsub)

        # calculate time-averaged velocity profiles
        UMeanAvg = np.dot( dtsub, self.U_mean[js:je,:] ) / np.sum(dtsub)
        VMeanAvg = np.dot( dtsub, self.V_mean[js:je,:] ) / np.sum(dtsub)
        WMeanAvg = np.dot( dtsub, self.W_mean[js:je,:] ) / np.sum(dtsub)

        # calculate time-averaged variances
        uuMeanAvg = np.dot( dtsub, self.uu_mean[js:je,:] ) / np.sum(dtsub)
        vvMeanAvg = np.dot( dtsub, self.vv_mean[js:je,:] ) / np.sum(dtsub)
        uvMeanAvg = np.dot( dtsub, self.uv_mean[js:je,:] ) / np.sum(dtsub)
        wwMeanAvg = np.dot( dtsub, self.ww_mean[js:je,:] ) / np.sum(dtsub)
        if SFS:
            if verbose: print 'Adding SFS component'
            uuMeanAvg += np.dot( dtsub, self.R11_mean[js:je,:] ) / np.sum(dtsub)
            vvMeanAvg += np.dot( dtsub, self.R22_mean[js:je,:] ) / np.sum(dtsub)
            uvMeanAvg += np.dot( dtsub, self.R12_mean[js:je,:] ) / np.sum(dtsub)
            wwMeanAvg += np.dot( dtsub, self.R33_mean[js:je,:] ) / np.sum(dtsub)

        # interpolate values to input heights
        Ux = np.interp( heights, self.hLevelsCell, UMeanAvg )
        Uy = np.interp( heights, self.hLevelsCell, VMeanAvg )
        Uz = np.interp( heights, self.hLevelsCell, WMeanAvg )
        sqrt_uuMeanAvg = np.interp( heights, self.hLevelsCell, np.sqrt(uuMeanAvg) ) # for agreement with variances_avg_cell.m
        #sqrt_uvMeanAvg = np.interp( heights, self.hLevelsCell, np.sqrt(uvMeanAvg) ) # not used
        sqrt_vvMeanAvg = np.interp( heights, self.hLevelsCell, np.sqrt(vvMeanAvg) )
        sqrt_wwMeanAvg = np.interp( heights, self.hLevelsCell, np.sqrt(wwMeanAvg) )
        uuMeanAvg = np.interp( heights, self.hLevelsCell, uuMeanAvg )
        uvMeanAvg = np.interp( heights, self.hLevelsCell, uvMeanAvg )
        vvMeanAvg = np.interp( heights, self.hLevelsCell, vvMeanAvg )
        wwMeanAvg = np.interp( heights, self.hLevelsCell, wwMeanAvg )

        Umag = np.sqrt( Ux**2 + Uy**2 + Uz**2 )
        if verbose:
            for zi,ui in zip(heights,Umag): print 'Umag at z=',zi,'m : ',ui,'m/s'

        # calculate wind direction
        windDir = np.abs( np.arctan2(Uy,Ux) )

        # calculate TKE and TI
        TIx = sqrt_uuMeanAvg / Umag #np.sqrt( uuMeanAvg ) / Umag
        TIy = sqrt_vvMeanAvg / Umag #np.sqrt( vvMeanAvg ) / Umag
        TIz = sqrt_wwMeanAvg / Umag #np.sqrt( wwMeanAvg ) / Umag
        TKE = 0.5*( uuMeanAvg + vvMeanAvg + wwMeanAvg )
        TIxyz = np.sqrt( 2./3.*TKE ) / Umag

        TIdir = uuMeanAvg *   np.cos(windDir)**2 \
              + uvMeanAvg * 2*np.sin(windDir)*np.cos(windDir) \
              + vvMeanAvg *   np.sin(windDir)**2
        TIdir = np.sqrt(TIdir) / Umag

        # save attributes
        self.TIx    = TIx
        self.TIy    = TIy
        self.TIz    = TIz
        self.TIdir  = TIdir
        self.TIxyz  = TIxyz
        self.TKE    = TKE

        if Nout==1:
            return TIx[0],TIy[0],TIz[0],TIdir[0],TIxyz[0],TKE[0]
        # }}}
        return TIx, TIy, TIz, TIdir, TIxyz, TKE

    def calcTI_hist(self,
            heights=[],
            tavg_window=600.0,
            dt=1.0,
            SFS=True,
            verbose=True):
        """ Calculate the turbulence intensity (TI) of the resolved
        fluctuations alone or combined fluctuations (including resolved
        and sub-filter scale, SFS).

        INPUTS
            tavg_window     size of window for moving average [s]
            heights         vertical locations at which to calculate TI [m]
            dt              uniform time interval to which to interpolate [s]
            SFS             set to True to include SFS terms

        OUTPUTS
            tavg            uniformly spaced times at which a moving average was calculated
            TIx_hist        variance in x-dir, TI*.shape = ( len(tavg), len(heights) )
            TIy_hist        variance in y-dir
            TIz_hist        variance in z-dir
            TIdir_hist      variance resolved to flow direction
            TIxyz_hist      variance assuming homogeneous turbulence, calculated from TKE
            TKE_hist        turbulent kinetic energy (TKE)
        """# {{{
        try:
            from scipy.ndimage import uniform_filter
        except ImportError:
            print 'Moving average calculation uses scipy.ndimage'
            return

        self.getVarsIfNeeded('uu_mean','vv_mean','ww_mean','uv_mean','uw_mean','vw_mean')
        if SFS: self.getVarsIfNeeded('R11_mean','R22_mean','R33_mean','R12_mean','R13_mean','R23_mean')

        Nout  = len(heights)
        if Nout==0:
            print 'Need to specify output heights'
            return

        # check for inconsistent array lengths
        Nt0 = len(self.t)
        fieldsToCheck = [
            self.U_mean, self.V_mean, self.W_mean,
            self.uu_mean, self.uv_mean, self.vv_mean, self.ww_mean,
            self.R11_mean, self.R12_mean, self.R22_mean, self.R33_mean
        ]
        if any([ field.shape[0] < Nt0 for field in fieldsToCheck ]):
            # need to prune arrays
            Nt_new = np.min([ field.shape[0] for field in fieldsToCheck ])
            print 'Inconsistent averaging field lengths... is simulation still running?'
            print '  truncated field histories from',Nt0,'to',Nt_new
            self.t = self.t[:Nt_new]
            self.U_mean = self.U_mean[:Nt_new,:]
            self.V_mean = self.V_mean[:Nt_new,:]
            self.W_mean = self.W_mean[:Nt_new,:]
            self.uu_mean = self.uu_mean[:Nt_new,:]
            self.uv_mean = self.uv_mean[:Nt_new,:]
            self.vv_mean = self.vv_mean[:Nt_new,:]
            self.ww_mean = self.ww_mean[:Nt_new,:]
            self.R11_mean = self.R11_mean[:Nt_new,:]
            self.R12_mean = self.R12_mean[:Nt_new,:]
            self.R22_mean = self.R22_mean[:Nt_new,:]
            self.R33_mean = self.R33_mean[:Nt_new,:]

        # setup uniform points for interpolation and averaging windows
        Nt = int(np.ceil((self.t[-1]-self.t[0])/dt))
        tuniform = np.arange(1,Nt+1)*dt + self.t[0]
        Navg    = int(tavg_window/dt)
        tavg    = tuniform[Navg/2:-Navg/2+1]
        Ntavg   = len(tavg)
        if verbose:
            print 'Interpolating to',Nt,'uniformly-spaced data points'
            print 'Moving average window:',tavg_window,'s'

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
            UMeanAvg = uniform_filter( U_mean_uniform, Navg )[Navg/2:-Navg/2+1] # length=Ntavg
            VMeanAvg = uniform_filter( V_mean_uniform, Navg )[Navg/2:-Navg/2+1]
            WMeanAvg = uniform_filter( W_mean_uniform, Navg )[Navg/2:-Navg/2+1]

            # calculate time-averaged variances
            uu_mean_interp = self.uu_mean[:,k-1] + frac*(self.uu_mean[:,k] - self.uu_mean[:,k-1]) # length=len(self.t)
            vv_mean_interp = self.vv_mean[:,k-1] + frac*(self.vv_mean[:,k] - self.vv_mean[:,k-1])
            uv_mean_interp = self.uv_mean[:,k-1] + frac*(self.uv_mean[:,k] - self.uv_mean[:,k-1])
            ww_mean_interp = self.ww_mean[:,k-1] + frac*(self.ww_mean[:,k] - self.ww_mean[:,k-1])
            uu_mean_uniform = np.interp( tuniform, self.t, uu_mean_interp ) # length=Nt
            vv_mean_uniform = np.interp( tuniform, self.t, vv_mean_interp )
            uv_mean_uniform = np.interp( tuniform, self.t, uv_mean_interp )
            ww_mean_uniform = np.interp( tuniform, self.t, ww_mean_interp )
            uuMeanAvg = uniform_filter( uu_mean_uniform, Navg )[Navg/2:-Navg/2+1] # length=Ntavg
            vvMeanAvg = uniform_filter( vv_mean_uniform, Navg )[Navg/2:-Navg/2+1]
            uvMeanAvg = uniform_filter( uv_mean_uniform, Navg )[Navg/2:-Navg/2+1]
            wwMeanAvg = uniform_filter( ww_mean_uniform, Navg )[Navg/2:-Navg/2+1]
            if SFS:
                if verbose: print 'Adding SFS component'
                R11_mean_interp = self.R11_mean[:,k-1] + frac*(self.R11_mean[:,k] - self.R11_mean[:,k-1]) # length=len(self.t)
                R22_mean_interp = self.R22_mean[:,k-1] + frac*(self.R22_mean[:,k] - self.R22_mean[:,k-1])
                R12_mean_interp = self.R12_mean[:,k-1] + frac*(self.R12_mean[:,k] - self.R12_mean[:,k-1])
                R33_mean_interp = self.R33_mean[:,k-1] + frac*(self.R33_mean[:,k] - self.R33_mean[:,k-1])
                R11_mean_uniform = np.interp( tuniform, self.t, R11_mean_interp ) #length=Nt
                R22_mean_uniform = np.interp( tuniform, self.t, R22_mean_interp )
                R12_mean_uniform = np.interp( tuniform, self.t, R12_mean_interp )
                R33_mean_uniform = np.interp( tuniform, self.t, R33_mean_interp )
                uuMeanAvg += uniform_filter( R11_mean_uniform, Navg )[Navg/2:-Navg/2+1] # length=Ntavg
                vvMeanAvg += uniform_filter( R22_mean_uniform, Navg )[Navg/2:-Navg/2+1]
                uvMeanAvg += uniform_filter( R12_mean_uniform, Navg )[Navg/2:-Navg/2+1]
                wwMeanAvg += uniform_filter( R33_mean_uniform, Navg )[Navg/2:-Navg/2+1]

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

        # save attributes
        self.tavg       = tavg
        self.TIx_hist   = TIx
        self.TIy_hist   = TIy
        self.TIz_hist   = TIz
        self.TIdir_hist = TIdir
        self.TIxyz_hist = TIxyz
        self.TKE_hist   = TKE
        # }}}
        return tavg, TIx, TIy, TIz, TIdir, TIxyz, TKE

    def calcShear(self,
            heights=[20.0,40.0,80.0],
            zref=80.0,
            Uref=8.0,
            interp=False,
            verbose=True):
        """ Estimate the shear from the average streamwise velocity from
        the final time step. Sets self.approxWindProfile to the fitted
        wind profile.

        INPUTS
            heights     list of three points to use to fit the wind profile

        OUTPUTS
            alpha       power law wind profile exponent
        """#{{{
        from scipy.interpolate import interp1d
        Uh = np.sqrt( self.U_mean**2 + self.V_mean**2 )[-1,:] # horizontal wind

        if interp:
            Ufn = interp1d( self.hLevelsCell, Uh, kind='linear' )
            U = Ufn(heights)
            self.approxHeights = heights
            self.approxU = U # at heights
            self.approxUfn = Ufn
        else:
            # find nearest cell without interpolating
            idx = [ np.argmin(np.abs(h-self.hLevelsCell)) for h in sorted(heights) ]
            assert(len(set(idx)) >= 3)
            heights = self.hLevelsCell[idx]
            U = Uh[idx]

        if verbose:
            print 'Estimating shear coefficient for Uref=',Uref,'and zref=',zref,':'
            print '     U=',U,'m/s'
            print '  at z=',heights,'m'
        lnz = np.log( np.array(heights)/zref )
        lnU = np.log( U/Uref )
        alpha = lnz.dot(lnU) / lnz.dot(lnz)

        self.Uh = Uh
        self.approxWindProfile = Uref * (self.hLevelsCell/zref)**alpha
        self.alpha = alpha

        #}}}
        return self.alpha

    def calcVeer(self, zhub=80.0, D=126.0, verbose=True):
        """ Estimate the veer from the average streamwise velocity from
        the final time step. Also calculates the height-varying wind
        direction, saved in self.windDir [deg].

        INPUTS
            zhub        hub height [m]
            D           rotor diameter [m]

        OUTPUTS
            veer        veer angle [deg]
                        >0 implies clockwise change in wind direction seen from above
        """#{{{
        dir = np.arctan2( -self.U_mean, -self.V_mean )[-1,:]
        dir[dir<0] += 2*np.pi
        dir *= 180.0/np.pi

        self.windDir = dir
        
        idxs = (self.hLevelsCell >= zhub-D/2) & (self.hLevelsCell <= zhub+D/2)
        if verbose:
            z = self.hLevelsCell[idxs]
            #for zi,angi,tf in zip(self.hLevelsCell,dir,idxs):
            #    print zi,'m ',angi,'deg ',tf
            print 'Estimating shear between z=',z[0],'and',z[-1],'m'
        rotorWindDir = dir[idxs]
        self.veer = rotorWindDir[-1] - rotorWindDir[0]

        #}}}
        return self.veer

    def calcTGradients(self, zi):
        """ Calculate the temperature gradient at the inversion height
        and at the top of the computational domain.

        INPUTS
            zi          inversion height [m]

        OUTPUTS
            dTdz_inv    temperature gradient at the inversion [deg K/m]
            dTdz_upper  upper temperature gradient [deg K/m]
        """#{{{
        self.getVarsIfNeeded('T_mean')

        ii = np.argmin(np.abs(self.hLevelsCell-zi))
        dTdz_inv = (self.T_mean[-1,ii+1] - self.T_mean[-1,ii]) \
                 / (self.hLevelsCell[ii+1] - self.hLevelsCell[ii])

        dTdz_upper = (self.T_mean[-1,-1] - self.T_mean[-1,-2]) \
                   / (self.hLevelsCell[-1] - self.hLevelsCell[-2])
#}}}
        return dTdz_inv, dTdz_upper

    def calcRichardsonNumber(self, g=9.81, zref=90.0, D=126.0, verbose=True):
        """ Estimate the Richardson number from the averaged T profile
        at the final time step which should range between -O(0.01) to
        +O(0.1), for unstable to stable. Difference formula used are
        second-order accurate.

        INPUTS
            g           gravity [m/s^2]
            zref        reference height used to locate the top of the rotor [m]
            D           rotor diameter used to locate the top of the rotor [m]

        OUTPUTS
            Ri          Richardson number
        """#{{{
        self.getVarsIfNeeded('T_mean')

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
            print 'Calculating Ri with:'
            print '  mean T :',Tmean
            print '  dT/dz at z=',z[1],':',dTdz1,' (finite difference)'
            print '  dT/dz at z=',z[0],':',dTdz0,' (finite difference)'
            print '  dT/dz at z=0:',dTdz,' (extrapolated)'
            print '  dU/dz at z=',z[1],':',dUdz1,' (finite difference)'
            print '  dU/dz at z=',z[0],':',dUdz0,' (finite difference)'
            print '  dU/dz at z=0:',dUdz,' (extrapolated)'

        self.Tsurf = Tsurf
        self.dUdz_surf = dUdz
        self.dTdz_surf = dTdz
        self.Ri = g/Tmean * dTdz / dUdz**2

        #}}}
        return self.Ri

#===========================================================
if __name__ == '__main__':
    import glob

    subdirs = glob.glob('*')
    data = read(*subdirs)

    TIx,TIy,TIz,TIdir,TIxyz,TKE = data.TI(avg_time=15000.,heights=[90.],SFS=True)
    print 'TIx=',TIx
    print 'TIy=',TIy
    print 'TIz=',TIz
    print 'TIdir=',TIdir
    print 'TIxyz=',TIxyz
    print 'TKE=',TKE

    tavg,TIx_hist,TIy_hist,TIz_hist,TIdir_hist,TIxyz_hist,TKE_hist = data.TI_hist(heights=[90.],SFS=True)

