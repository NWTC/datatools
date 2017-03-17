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

def readRotorOutputs(datadir='turbineOutput',Nturb=1,turbineList=None):
    """Returns a dictionary of pandas dataframes for each turbine

    If turbineList is specified, then selected turbines are returned; 
    otherwise, all turbines (presumably Nturb many) are returned.
    """
    dataseries = TimeSeries(datadir,verbose=False)
    rotorDataNames = dataseries.outputs('rotor')
    turbinedata = dict()

    if turbineList is None:
        turbineList = range(Nturb)

    for dataname in rotorDataNames:
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
                # TODO: fix header names for data files with blade points
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

