import numpy as np
import pandas as pd

from SOWFA.timeseries import TimeSeries

def _processTurbineOutputHeader(line):
    """turbineOutput file headers have the following format:

    #Turbine    [Blade]    Time(s)    dt(s)    outputQuantity (units)

    where the 'Blade' column only occurs for blade-sampled output. 
    In this case, 'outputQuantity' has 'numBladePoints' (defined in
    turbineArrayProperties) points.
    """
    if line.startswith('#'):
        line = line[1:]
    line = line.split()
    headerNames = []
    for i,name in enumerate(line):
        headerNames.append(name)
        if name.startswith('dt'):
            break
    headerNames.append(' '.join(line[i+1:]))
    return headerNames

def readRotorOutputs(datadir='turbineOutput',prefix='rotor',Nturb=1,turbineList=None):
    """Returns a dictionary of pandas dataframes for each turbine

    If turbineList is specified, then selected turbines are returned; 
    otherwise, all turbines (presumably Nturb many) are returned.
    """
    dataseries = TimeSeries(datadir,verbose=False)
    #dataNames = dataseries.outputs(prefix)
    dataNames = dataseries.outputs(prefix) + ['bladePitch']
    turbinedata = dict()

    if turbineList is None:
        turbineList = range(Nturb)

    for dataname in dataNames:
        print 'Processing',dataname
        dataseries.setFilename(dataname)
        
        dframes = []
        
        # loop over restarts
        for irest,fname in enumerate(dataseries):
            print '  datafile',irest,':',fname
            with open(fname,'r') as f:
                headerNames = _processTurbineOutputHeader(f.readline())
                next_df = pd.read_csv(f,
                                      delim_whitespace=True,
                                      header=None, names=headerNames, 
                                     ).set_index('Time(s)') # Note: setting the index removes the time column
            
            for iturb,turbNum in enumerate(turbineList):
                next_df = next_df.loc[next_df['Turbine'] == turbNum]
                if irest==0:
                    dframes.append(next_df)
                    assert(len(dframes) == iturb+1)
                else:
                    dframes[iturb] = next_df.combine_first(dframes[iturb]) # if overlap, overwrite with next_df values
        
        # append full time series (with all restarts) to the complete turbinedata frame for each turbine
        for iturb,turbNum in enumerate(turbineList):
            if turbNum not in turbinedata.keys():
                turbinedata[turbNum] = dframes[iturb]
            else:
                turbinedata[turbNum] = pd.concat(
                        (turbinedata[turbNum], dframes[iturb].iloc[:,len(headerNames)-2:]),
                        axis=1)

    # sort everything by the (time) index once
    for turbNum in turbineList:
        turbinedata[turbNum].sort_index(inplace=True)

    return turbinedata

def readBladeOutputs(datadir='turbineOutput',prefix='blade',Nturb=1,bladeNum=0,turbineList=None):
    """Returns a dictionary of pandas dataframes for each turbine

    If turbineList is specified, then selected turbines are returned; 
    otherwise, all turbines (presumably Nturb many) are returned.

    For reference: http://pandas.pydata.org/pandas-docs/stable/reshaping.html
    """
    dataseries = TimeSeries(datadir,verbose=False)
    dataNames = dataseries.outputs(prefix)
    dataNames.remove('bladePitch') # this is the collective pitch for all blades, and should be a 'rotor' quantity
    turbinedata = dict()

    if turbineList is None:
        turbineList = range(Nturb)

    for dataname in dataNames:
        print 'Processing',dataname
        dataseries.setFilename(dataname)
        
        dframes = []
        
        # loop over restarts
        for irest,fname in enumerate(dataseries):
            print '  datafile',irest,':',fname
            with open(fname,'r') as f:
                headerNames = _processTurbineOutputHeader(f.readline())
                testline = f.readline().split()
                numBladePoints = len(testline) - len(headerNames) + 1
                print '  (detected',numBladePoints,'blade points)'
                fieldname = headerNames[-1]
                headerNames = headerNames[:-1] + [ fieldname+'_'+str(ipt) for ipt in range(numBladePoints) ]
            with open(fname,'r') as f:
                f.readline() # skip header
                next_df = pd.read_csv(f,
                                      delim_whitespace=True,
                                      header=None, names=headerNames, 
                                     ).set_index('Time(s)') # Note: setting the index removes the time column
            
            for iturb,turbNum in enumerate(turbineList):
                next_df = next_df.loc[(next_df['Turbine'] == turbNum) & (next_df['Blade'] == bladeNum)]
                if irest==0:
                    dframes.append(next_df)
                    assert(len(dframes) == iturb+1)
                else:
                    dframes[iturb] = next_df.combine_first(dframes[iturb]) # if overlap, overwrite with next_df values
        
        # append full time series (with all restarts) to the complete turbinedata frame for each turbine
        for iturb,turbNum in enumerate(turbineList):
            if turbNum not in turbinedata.keys():
                turbinedata[turbNum] = dframes[iturb]
            else:
                turbinedata[turbNum] = pd.concat(
                        (turbinedata[turbNum], dframes[iturb].iloc[:,len(headerNames)-2:]),
                        axis=1)

    # sort everything by the (time) index once
    for turbNum in turbineList:
        turbinedata[turbNum].sort_index(inplace=True)

    return turbinedata

