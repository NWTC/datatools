"""
WFIP 2 Data Processing Tools
============================
Written by Eliot Quon (eliot.quon@nrel.gov)

Tools to read series of data files contained in either a single
directory, or a series of subdirectories, into a single pandas
dataframe.

Also includes helper tools, e.g., to make time-height wind plots.

Sample usage:
    from datatools import wfip2, remote_sensing
    # to read $dpath/*.winds
    df = wfip2.read_dir(dpath,
                        reader=remote_sensing.ESRL_wind_profiler,
                        ext='winds',
                        na_values=[999999,-980.0])

"""
from __future__ import print_function
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    # use NOAA colormap
    from abl.miscellaneous import idl_colortable
except ImportError:
    windspeed_colormap = plt.cm.gist_ncar
else:
    windspeed_colormap = idl_colortable()

winddirection_colormap = plt.cm.hsv

#
# Wrappers for data loading
#

def read_dir(dpath='.',
             reader=pd.read_csv,
             file_filter='*',
             ext='csv',
             sort=True,
             verbose=False,
             **kwargs):
    """Wrapper around pandas read_csv() or a custom remote_sensing.*
    reader function
    
    Returns concatenated dataframe made up of dataframes read from CSV
    files in specified directory

    Additional keyword arguments are passed to the data reader.
    """
    dataframes = []
    fpathlist = glob.glob(os.path.join(dpath,file_filter))
    if sort:
        fpathlist.sort()
    for fpath in fpathlist:
        if not os.path.isfile(fpath): continue
        if not fpath.endswith(ext):
            continue
        #fname = os.path.split(fpath)[-1]
        #if verbose: print('Reading '+fname)
        if verbose: print('Reading '+fpath)
        df = reader(fpath,**kwargs)
        dataframes.append(df)
    return pd.concat(dataframes)


def read_date_dirs(dpath='.',
                   reader=pd.read_csv,
                   ext='csv',
                   expected_date_format=None,
                   verbose=False,
                   **kwargs):
    """Wrapper around pandas read_csv() or custom remote_sensing.*
    reader function
    
    Return concatenated dataframe made up of dataframes read from
    CSV files contained in date subdirectories. 

    Additional keyword arguments are passed to the data reader.
    """
    dataframes = []
    for dname in os.listdir(dpath):
        Nfiles = 0
        fullpath = os.path.join(dpath,dname)
        if os.path.isdir(fullpath):
            try:
                collection_date = pd.to_datetime(dname,
                                                 format=expected_date_format)
            except ValueError:
                if verbose: print('Skipping '+dname)
            else:
                print('Processing '+fullpath)
                for fname in os.listdir(fullpath):
                    fpath = os.path.join(fullpath,fname)
                    if not fname.endswith(ext): continue
                    if verbose: print('  reading '+fname)
                    df = reader(fpath,**kwargs)
                    dataframes.append(df)
                    Nfiles += 1
            print('  {} dataframes added'.format(Nfiles))
    return pd.concat(dataframes)


#
# Time-height plotting tools
#

def plot_wind(df,
              height_name='height',
              speed_name='speed',
              direction_name='direction',
              datetime_range=(None,None),
              verbose=False):
    """Make time-height plot of wind speed and direction, assuming a
    datetime index has been set
    """
    # set up time range
    if datetime_range[0] is None:
        tstart = df.index[0]
    else:
        tstart = pd.to_datetime(datetime_range[0])
    if datetime_range[1] is None:
        tend = df.index[-1]
    else:
        tend = pd.to_datetime(datetime_range[1])
    trange = (df.index >= tstart) & (df.index <= tend)

    # get wind history subset
    dfsub = df.loc[trange]
    height = dfsub[height_name].unique()
    time = dfsub.loc[dfsub[height_name]==height[0]].index
    if verbose:
        print('heights:',height)
        print('times:',time)
    X,Y = np.meshgrid(time.to_pydatetime(),height,indexing='ij')
    Nt, Nh = len(time), len(height)
    wspd = np.zeros((Nt,Nh))
    wdir = np.zeros((Nt,Nh))
    for k,h in enumerate(height):
        wspd[:,k] = dfsub.loc[dfsub[height_name]==h,speed_name]
        wdir[:,k] = dfsub.loc[dfsub[height_name]==h,direction_name]

    # make plot
    fig,ax = plt.subplots(nrows=2,sharex=True,sharey=True,figsize=(10,6))
    wslevels = np.arange(0,25.1,0.5)
    wdlevels = np.arange(0,360.1,7.5)

    cont = ax[0].contourf(X, Y, wspd, levels=wslevels, cmap=windspeed_colormap)
    cbar = fig.colorbar(cont, ax=ax[0], ticks=np.arange(0,26),
                        label='wind speed [m/s]')

    cont = ax[1].contourf(X, Y, wdir, levels=wdlevels, cmap=winddirection_colormap)
    cbar = fig.colorbar(cont, ax=ax[1], ticks=np.arange(0,361,45),
                        label='wind direction [deg]')

    return fig, ax

def plot_temp(df,
              height_name='height',
              temperature_name='temperature',
              datetime_range=(None,None),
              contour_res=0.5,
              verbose=False):
    """Make time-height plot of temperature, assuming a datetime index
    has been set
    """
    # set up time range
    if datetime_range[0] is None:
        tstart = df.index[0]
    else:
        tstart = pd.to_datetime(datetime_range[0])
    if datetime_range[1] is None:
        tend = df.index[-1]
    else:
        tend = pd.to_datetime(datetime_range[1])
    trange = (df.index >= tstart) & (df.index <= tend)

    # get wind history subset
    dfsub = df.loc[trange]
    height = dfsub[height_name].unique()
    time = dfsub.loc[dfsub[height_name]==height[0]].index
    if verbose:
        print('heights:',height)
        print('times:',time)
    X,Y = np.meshgrid(time.to_pydatetime(),height,indexing='ij')
    Nt, Nh = len(time), len(height)
    thetav = np.zeros((Nt,Nh))
    for k,h in enumerate(height):
        thetav[:,k] = dfsub.loc[dfsub[height_name]==h,temperature_name]

    # make plot
    fig,ax = plt.subplots(sharex=True,sharey=True,figsize=(10,3))
    tlevels = np.arange(np.round(np.nanmin(thetav)/contour_res)*contour_res,
                        np.round(np.nanmax(thetav)/contour_res)*contour_res+0.1,
                        contour_res)
    cont = ax.contourf(X,Y,thetav, levels=tlevels, cmap='inferno')
    cbar = fig.colorbar(cont, label='potential temperature [K]')

    return fig, ax


#
# Profile extraction
#

def get_profile_at_time(df,time,field='speed',height_name='height'):
    """Interpolate field data to specified time, assuming a datetime
    index has been set.

    Returns height and data vectors
    """
    wide = df.pivot(columns=height_name,values=field)
    wide.loc[time] = None
    wide = wide.interpolate(method='slinear')
    profile = wide.loc[time]
    return profile.index, profile.values

def average_profile_over_times(df,trange,field='speed',height_name='height'):
    """Temporal average of field data over specified time range,
    assuming a datetime index has been set.

    Returns height and data vectors
    """
    trange = (df.index >= trange[0]) & (df.index <= trange[1])
    dfsub = df.loc[trange]
    times = dfsub.index.unique()
    print('Average over',len(times),times)
    wide = dfsub.pivot(columns=height_name,values=field)
    profile = wide.mean(axis=0)
    z = dfsub[height_name].unique()
    return z, profile.values

def stdev_profile_over_times(df,trange,field='speed',height_name='height'):
    """Standard deviation of field data over specified time range,
    assuming a datetime index has been set.

    Returns height and data vectors
    """
    trange = (df.index >= trange[0]) & (df.index <= trange[1])
    dfsub = df.loc[trange]
    times = dfsub.index.unique()
    print('Standard deviation over',len(times),times)
    wide = dfsub.pivot(columns=height_name,values=field)
    profile = wide.std(axis=0)
    z = dfsub[height_name].unique()
    return z, profile.values

