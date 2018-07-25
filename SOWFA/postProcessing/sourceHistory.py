"""
For processing SOWFA driving force data in postProcessing/SourceHistory
Based on averaging.py written by Eliot Quon

written by Dries Allaerts (dries.allaerts@nrel.gov)

The class can handle both height-dependent and constant source files

Sample usage:

    from datatools.SOWFA.postProcessing.sourceHistory import SourceHistory

    # read all time directories in current working directory or in subdirectory called 'SourceHistory'
    srcData = SourceHistory()

    # read '0' and '1000' in current working directory
    srcData = SourceHistory( 0, 1000 )

    # read all time directories in specified directory
    srcData = SourceHistory('caseX/postProcessing/SourceHistory')

    # read specified time directories
    srcData = SourceHistory('caseX/postProcessing/SourceHistory/0',
                        'caseX/postProcessing/SourceHistory/1000')

    # read UX source term in all time directories in current working directory
    srcData = SourceHistory(varList='UX')

"""
from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt

sources = ['UX','UY','UZ','T']

class SourceHistory(object):

    def __init__(self,*args,**kwargs):
        """Find and process all time directories"""
        self._processed = []
        self.hLevelsCell = None
        self.simTimeDirs = [] # output time names
        self.simStartTimes = [] # start or restart simulation times
        self.imax = None # for truncating time series
        
        if len(args)==0:
            args = os.listdir('.')
            if 'SourceHistory' in args: args = ['SourceHistory']

        # find results
        for arg in args:
            if isinstance(arg,str):
                if not os.path.isdir(arg): continue
            else: #float or int were given
                arg = '{:g}'.format(arg)
                if not os.path.isdir(arg): continue

            if arg[-1] == os.sep: arg = arg[:-1] # strip trailing slash

            listing = os.listdir(arg)
            if 'SourceTHistory' in listing:
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
                        if 'SourceTHistory' in os.listdir(arg+os.sep+dirname):
                            self.simTimeDirs.append( arg+os.sep+dirname )
                            self.simStartTimes.append( startTime )
                    except ValueError: pass # dirname is not a number

        # sort results
        self.simTimeDirs = [ x[1] for x in sorted(zip(self.simStartTimes,self.simTimeDirs)) ]
        self.simStartTimes.sort()

        print('Simulation (re)start times:',self.simStartTimes)

        # process all output dirs
        if len(self.simTimeDirs) > 0:
            self._processdirs( self.simTimeDirs, **kwargs )
        else:
            print('No averaging time directories found!')
    

    def __repr__(self):
        s = 'SOWFA postProcessing: driving force data'
        for t,d in zip(self.simStartTimes,self.simTimeDirs):
            fullpath = os.path.realpath(os.curdir) + os.sep + d
            s += '\n  {:f}\t{:s}'.format(t,fullpath)
        s += '\ntimes: {:d} [{:f},{:f}]'.format(len(self.t),self.t[0],self.t[-1])
        s += '\nheights: {:d} [{:f},{:f}]'.format(len(self.hLevelsCell),self.hLevelsCell[0],self.hLevelsCell[-1])
        return s

    def _processdirs(self,tdirList,varList=['UX','UY','UZ','T']):
        """Reads all files within an output time directory.
        An object attribute corresponding to a source term
        output name is updated, e.g.:
            ${timeDir}/SourceUXHistory is appended to the array self.UX

        Typically, objects have shape (Nt,Nz).
        """
        outputs = []
        if isinstance( varList, (str,) ):
            if varList.lower()=='all':
                # special case: read all vars
                allOutputs = os.listdir(tdirList[0])
                for field in allOutputs:
                    outputs.append( field.lstrip('Source').rstrip('History') )
            else: # specified single var
                outputs = [varList]
        else: # specified list
            outputs = varList

        # Check whether sources depend on height or not
        with open(tdirList[0]+os.sep+'Source'+outputs[0]+'History','r') as f:
            line = f.readline().split()
        if line[0] == 'Time':
            self.hLevelsCell = np.array([0.0])
            header_lines = 1
        elif line[0] == 'Heights':
            self.hLevelsCell = np.array([ float(val) for val in line[2:] ])
            header_lines = 2
        else:
            print('Error: Expected first line to start with "Time" or "Heights", but instead read',line[0])
            return

        # check that we have the same amount of data
        for tdir in tdirList:
            Nlines = []
            for field in outputs:
                output = tdir + os.sep + 'Source'+field+'History'
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
        for field in outputs:
            arrays = [ np.genfromtxt( tdir+os.sep+'Source'+field+'History', skip_header = header_lines) for tdir in tdirList ]
            newdata = np.concatenate(arrays)
            setattr( self, field, newdata[:self.imax,2:] )

            self._processed.append(field)
            print('  read',field)

        self.t = np.array( newdata[:self.imax,0] )
        self.dt = np.array( newdata[:self.imax,1] )
        
        return None

    def _trim_series_if_needed(self,fields_to_check=sources):
        """check for inconsistent array lengths and trim if needed"""
        Nt0 = len(self.t)
        for field in fields_to_check:
            try:
                getattr(self,field)
            except AttributeError:
                fields_to_check.remove(field)
        Nt_new = np.min([ getattr(self,field).shape[0] for field in fields_to_check ])
        if Nt_new < Nt0:
            # need to prune arrays
            print('Inconsistent averaging field lengths... is simulation still running?')
            print('  truncated field histories from',Nt0,'to',Nt_new)
            self.t = self.t[:Nt_new]
            for field in fields_to_check:
                setattr(self, field, getattr(self,field)[:Nt_new,:])
            self.imax = Nt_new
