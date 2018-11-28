#
# Tools for processing SOWFA actuator line output data located in
#   postProcessing/turbineOutput
#
# USAGE:
#   tdata = TurbineOutput()
#   tdata.
#
from __future__ import print_function
import numpy as np
import pandas as pd

from datatools.series import SOWFATimeSeries

class TurbineOutput(object):
    """Container for SOWFA turbineOutput data"""

    def __init__(self,datadir='turbineOutput',
                 Nturb=1,turbineList=None,
                 timeIndexName='Time(s)',toffset=None):
        """If turbineList is specified, then selected turbines are
        returned; otherwise, all turbines (presumably Nturb many) are
        returned.
        """
        self.datadir = datadir
        self.Nturb = Nturb
        if turbineList is None:
            turbineList = range(Nturb)
        self.turbineList = turbineList
        self.timeIndexName = timeIndexName
        self.toffset = toffset

    def processTurbineOutputHeader(self,line):
        """turbineOutput file headers have the following format:

        #Turbine    [Blade]    Time(s)    dt(s)    outputQuantity (units)

        where the 'Blade' column only occurs for blade-sampled output. 
        In this case, 'outputQuantity' has 'numBladePoints' (defined in
        turbineArrayProperties) points.
        """
        if line.startswith('#'):
            line = line.lstrip('#')
        line = line.split()
        headerNames = []
        for i,name in enumerate(line):
            headerNames.append(name)
            if name.startswith('dt'):
                break
        headerNames.append(' '.join(line[i+1:]))
        return headerNames

    def readRotorOutputs(self,prefix='rotor'):
        """Returns a dictionary of pandas dataframes for each turbine"""
        dataseries = SOWFATimeSeries(self.datadir,verbose=False)
        #dataNames = dataseries.outputs(prefix)
        dataNames = dataseries.outputs(prefix) + ['bladePitch']
        turbinedata = dict()

        for dataname in dataNames:
            print('Processing',dataname)
            dataseries.get(dataname)
            
            dframes = []
            
            # loop over restarts
            for irest,fname in enumerate(dataseries):
                print('  datafile',irest,':',fname)
                with open(fname,'r') as f:
                    headerNames = processTurbineOutputHeader(f.readline())
                    df = pd.read_csv(f,
                            delim_whitespace=True,
                            header=None, names=headerNames)

                if self.toffset is not None:
                    df[self.timeIndexName] -= self.toffset
                if self.timeIndexName is not None:
                    # Note: setting the index removes the time column
                    df.set_index(self.timeIndexName,inplace=True)

                for iturb,turbNum in enumerate(self.turbineList):
                    next_df = df.loc[df['Turbine'] == turbNum]
                    if irest==0:
                        dframes.append(next_df)
                        assert(len(dframes) == iturb+1)
                    else:
                        # if overlap, overwrite with next_df values
                        dframes[iturb] = next_df.combine_first(dframes[iturb])
            
            # append full time series (with all restarts) to the complete
            # turbinedata frame for each turbine
            for iturb,turbNum in enumerate(self.turbineList):
                if turbNum not in turbinedata.keys():
                    turbinedata[turbNum] = dframes[iturb]
                else:
                    turbinedata[turbNum] = pd.concat(
                            (turbinedata[turbNum], dframes[iturb].iloc[:,-1:]),
                            axis=1)

        # sort everything by the (time) index once
        for turbNum in self.turbineList:
            turbinedata[turbNum].sort_index(inplace=True)

        return turbinedata

    def readBladeOutputs(self,prefix='blade',bladeNum=0):
        """Returns a dictionary of pandas dataframes for each turbine

        For reference: http://pandas.pydata.org/pandas-docs/stable/reshaping.html
        """
        dataseries = SOWFATimeSeries(self.datadir,verbose=False)
        dataNames = dataseries.outputs(prefix)
        # Note: 'bladePitch' is the collective pitch for all blades, and should
        #       be a 'rotor' quantity
        dataNames.remove('bladePitch')
        turbinedata = dict()

        for dataname in dataNames:
            print('Processing',dataname)
            dataseries.get(dataname)
            
            dframes = []
            
            # loop over restarts
            for irest,fname in enumerate(dataseries):
                print('  datafile',irest,':',fname)
                with open(fname,'r') as f:
                    headerNames = processTurbineOutputHeader(f.readline())
                    testline = f.readline().split()
                    numBladePoints = len(testline) - len(headerNames) + 1
                    print('  (detected',numBladePoints,'blade points)')
                    fieldname = headerNames[-1]
                    headerNames = headerNames[:-1] \
                            + [ fieldname+'_'+str(ipt) for ipt in range(numBladePoints) ]
                with open(fname,'r') as f:
                    f.readline() # skip header
                    df = pd.read_csv(f,
                            delim_whitespace=True,
                            header=None, names=headerNames)
                
                if self.toffset is not None:
                    df[self.timeIndexName] -= self.toffset
                if self.timeIndexName is not None:
                    # Note: setting the index removes the time column
                    df.set_index(self.timeIndexName,inplace=True)

                for iturb,turbNum in enumerate(self.turbineList):
                    sel = (df['Turbine'] == turbNum) & (df['Blade'] == bladeNum)
                    next_df = df.loc[sel]
                    if irest==0:
                        dframes.append(next_df)
                        assert(len(dframes) == iturb+1)
                    else:
                        # if overlap, overwrite with next_df values
                        dframes[iturb] = next_df.combine_first(dframes[iturb])
            
            # append full time series (with all restarts) to the complete
            # turbinedata frame for each turbine
            for iturb,turbNum in enumerate(self.turbineList):
                if turbNum not in turbinedata.keys():
                    turbinedata[turbNum] = dframes[iturb]
                else:
                    turbinedata[turbNum] = pd.concat(
                            (turbinedata[turbNum], dframes[iturb].iloc[:,-numBladePoints:]),
                            axis=1)

        # sort everything by the (time) index once
        for turbNum in self.turbineList:
            turbinedata[turbNum].sort_index(inplace=True)

        return turbinedata

    def readTowerOutputs(self,prefix='tower'):
        """Returns a dictionary of pandas dataframes for each turbine

        For reference: http://pandas.pydata.org/pandas-docs/stable/reshaping.html
        """
        dataseries = SOWFATimeSeries(self.datadir,verbose=False)
        dataNames = dataseries.outputs(prefix)
        turbinedata = dict()

        for dataname in dataNames:
            print('Processing',dataname)
            dataseries.get(dataname)
            
            dframes = []
            
            # loop over restarts
            for irest,fname in enumerate(dataseries):
                print('  datafile',irest,':',fname)
                with open(fname,'r') as f:
                    headerNames = processTurbineOutputHeader(f.readline())
                    testline = f.readline().split()
                    numTowerPoints = len(testline) - len(headerNames) + 1
                    print('  (detected',numTowerPoints,'tower points)')
                    fieldname = headerNames[-1]
                    headerNames = headerNames[:-1] \
                            + [ fieldname+'_'+str(ipt) for ipt in range(numTowerPoints) ]
                with open(fname,'r') as f:
                    f.readline() # skip header
                    df = pd.read_csv(f,
                            delim_whitespace=True,
                            header=None, names=headerNames)
                
                if self.toffset is not None:
                    df[self.timeIndexName] -= self.toffset
                if self.timeIndexName is not None:
                    # Note: setting the index removes the time column
                    df.set_index(self.timeIndexName,inplace=True)

                for iturb,turbNum in enumerate(self.turbineList):
                    next_df = df.loc[df['Turbine'] == turbNum]
                    if irest==0:
                        dframes.append(next_df)
                        assert(len(dframes) == iturb+1)
                    else:
                        # if overlap, overwrite with next_df values
                        dframes[iturb] = next_df.combine_first(dframes[iturb])
            
            # append full time series (with all restarts) to the complete
            # turbinedata frame for each turbine
            for iturb,turbNum in enumerate(self.turbineList):
                if turbNum not in turbinedata.keys():
                    turbinedata[turbNum] = dframes[iturb]
                else:
                    turbinedata[turbNum] = pd.concat(
                            (turbinedata[turbNum], dframes[iturb].iloc[:,-numTowerPoints:]),
                            axis=1)

        # sort everything by the (time) index once
        for turbNum in self.turbineList:
            turbinedata[turbNum].sort_index(inplace=True)

        return turbinedata

