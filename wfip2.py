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
import matplotlib.dates as mdates

try:
    # use NOAA colormap
    from abl.miscellaneous import idl_colortable
except ImportError:
    windspeed_colormap = plt.cm.gist_ncar
else:
    windspeed_colormap = idl_colortable()

try:
    # matplotlib version >= 3.0 
    cyclic_colormap = plt.cm.twilight
except AttributeError:
    cyclic_colormap = plt.cm.hsv

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
    for dname in sorted(os.listdir(dpath)):
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

    DEPRECATED--use plot(), or plot_windspeed() and plot_winddirection()
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

    cont = ax[1].contourf(X, Y, wdir, levels=wdlevels, cmap=cyclic_colormap)
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

    DEPRECATED--use plot(), or plot_temperature()
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

def _trim_datetime(df, datetime_range):
    if datetime_range == (None,None):
        return df
    # set up time range
    if datetime_range[0] is None:
        tstart = df.index[0]
    else:
        tstart = pd.to_datetime(datetime_range[0])
    if datetime_range[1] is None:
        tend = df.index[-1]
    else:
        tend = pd.to_datetime(datetime_range[1])
    return df.loc[(df.index >= tstart) & (df.index <= tend)]

def plot_timeheight(df, column, ax=None,
                    datetime_range=(None,None),
                    datetime_name='datetime',
                    height_name='height',
                    cmap='viridis',
                    label='',
                    cbar_ticks=None,
                    **kwargs):
    """Make time-height plot

    If the dataframe does not have a datetime index, then columns
    'datetime_name' and 'height_name' are used to form lists of unique
    datetimes and heights to plot.

    If the dataframe has a datetime index, then only 'height_name' is
    used to form a list of unique heights to plot.

    If the dataframe has a multiindex, then the first should be a
    datetime index and the second should be the height index.
    """
    # setup index
    if isinstance(df.index, pd.core.index.MultiIndex):
        # multindex
        height = df.index.levels[1]
        df = df.reset_index(level=1)
    elif isinstance(df.index, pd.DatetimeIndex):
        # single datetime index
        height = df[height_name].unique()
        df = df.reset_index()
    else:
        # default range index
        height = df[height_name].unique()
        df = df.set_index(datetime_name)

    # trim rows outside of range
    df = _trim_datetime(df,datetime_range)
    time = df.index.unique().to_pydatetime()

    # create arrays for plotting
    X,Y = np.meshgrid(time, height, indexing='ij')
    #F = np.zeros((len(time),len(height)))
    #for k,h in enumerate(height):
    #    F[:,k] = df.loc[df.index.levels[1]==h, column]
    F = df.pivot(columns=height_name, values=column).values

    # make plot
    if ax is None:
        _,ax = plt.subplots(figsize=(10,3))
    cm = ax.pcolormesh(X,Y,F,cmap=cmap,**kwargs)
    cbar = plt.colorbar(cm,ax=ax,label=label)
    if cbar_ticks is not None:
        cbar.set_ticks(cbar_ticks)

    # format time axis
    ax.set_xlabel('')
    ax.set_ylabel('height [m]')
    #ax.xaxis.set_major_locator(mdates.DayLocator())
    #ax.xaxis.set_minor_locator(mdates.HourLocator())
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('\n%Y %b %d'))
    #ax.xaxis.set_minor_formatter(mdates.DateFormatter('%HZ'))


def plot_windspeed(df, name='windspeed', components=None,
                   cmap='viridis', vmin=None, vmax=None,
                   label='wind speed [m/s]',
                   **kwargs):
    """Plot wind speed by calling plot_timeheight()
    
    'name' corresponds to a column within the dataframe. Otherwise,
    'components', a list of u,v[,w] velocity components, is used to
    calculate the wind speed with 'name'.

    See plot_timeheight() for general keyword arguments.
    """
    if name not in df.columns:
        if components is None:
            raise KeyError("Column '{:s}' not found; specify correct column name or list of velocity components".format(name))
        else:
            print('Calculating wind speed from {:s}'.format(str(components)))
            df[name] = np.sqrt((df[components]**2).sum(axis=1))
    elif components is not None:
        print("Components {:s} ignored because column '{:s}' exists".format(
                str(components), name))

    plot_timeheight(df, column=name,
                    cmap=cmap,label=label,vmin=vmin,vmax=vmax,
                    **kwargs)

def plot_winddirection(df, name='winddirection', components=None,
                       cmap='viridis', vmin=None, vmax=None,
                       label='wind direction [deg]',
                       cbar_ticks=None,
                       full_360=False,
                       **kwargs):
    """Plot wind direction by calling plot_timeheight()

    If 'full_360' is True, then a cyclic colormap is used and vmin/max
    are ignored.
    
    'name' corresponds to a column within the dataframe. Otherwise,
    'components', a list of u,v velocity components, is used to
    calculate the wind direction with 'name'.

    See plot_timeheight() for general keyword arguments.
    """
    if name not in df.columns:
        if components is None:
            raise KeyError("Column '{:s}' not found; specify correct column name or list of velocity components".format(name))
        else:
            print('Calculating wind direction from {:s}'.format(str(components)))
            df[name] = 180./np.pi * np.arctan2(-df[components[0]], -df[components[1]])
            df.loc[df[name] < 0, name] += 360
    elif components is not None:
        print("Components {:s} ignored because column '{:s}' exists".format(
                str(components), name))

    if full_360:
        plot_timeheight(df, column=name,
                        cmap=cyclic_colormap,vmin=0,vmax=360,
                        label=label,
                        cbar_ticks=np.arange(0,361,45),
                        **kwargs)
    else:
        plot_timeheight(df, column=name,
                        cmap=cmap,vmin=vmin,vmax=vmax,
                        label=label,
                        cbar_ticks=cbar_ticks,
                        **kwargs)


def plot_temperature(df, name='T',
                     cmap='inferno', vmin=None, vmax=None,
                     label='temperature [K]',
                     **kwargs):
    """Plot temperature by calling plot_timeheight()
    
    'name' corresponds to a column within the dataframe.

    See plot_timeheight() for general keyword arguments.
    """
    if name not in df.columns:
        raise KeyError("Column '{:s}' not found; specify correct column name".format(name))

    plot_timeheight(df, column=name,
                    cmap=cmap,label=label,vmin=vmin,vmax=vmax,
                    **kwargs)

def plot_all(df, ax=None,
             windspeed_name='windspeed',
             windspeed_label='wind speed [m/s]',
             winddirection_name='winddirection',
             winddirection_label='wind direction [deg]',
             temperature_name='T',
             temperature_label='temperature [K]',
             **kwargs):
    if ax is None:
        fig,ax = plt.subplots(nrows=3,sharex=True,sharey=True,figsize=(10,8))
    plot_windspeed(df,
                ax=ax[0], 
                name=windspeed_name,
                label=windspeed_label,
                **kwargs)
    plot_winddirection(df,
                ax=ax[1], 
                name=winddirection_name,
                label=winddirection_label,
                **kwargs)
    plot_temperature(df,
                 ax=ax[2], 
                 name=temperature_name,
                 label=temperature_label,
                 **kwargs)
    return fig,ax

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

def average_profile_over_times(df,trange,field='speed',height_name='height',
                               verbose=True):
    """Temporal average of field data over specified time range,
    assuming a datetime index has been set.

    Returns height and data vectors
    """
    trange = (df.index >= trange[0]) & (df.index <= trange[1])
    dfsub = df.loc[trange]
    times = dfsub.index.unique()
    if verbose: print('Average over',len(times),times)
    wide = dfsub.pivot(columns=height_name,values=field)
    profile = wide.mean(axis=0)
    z = dfsub[height_name].unique()
    return z, profile.values

def stdev_profile_over_times(df,trange,field='speed',height_name='height',
                             verbose=True):
    """Standard deviation of field data over specified time range,
    assuming a datetime index has been set.

    Returns height and data vectors
    """
    trange = (df.index >= trange[0]) & (df.index <= trange[1])
    dfsub = df.loc[trange]
    times = dfsub.index.unique()
    if verbose: print('Standard deviation over',len(times),times)
    wide = dfsub.pivot(columns=height_name,values=field)
    profile = wide.std(axis=0)
    z = dfsub[height_name].unique()
    return z, profile.values

