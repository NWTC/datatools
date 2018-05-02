from __future__ import print_function
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    # use NOAA colormap
    from abl.meteorology import idl_colortable
except ImportError:
    windspeed_colormap = plt.cm.gist_ncar
else:
    windspeed_colormap = idl_colortable()

winddirection_colormap = plt.cm.hsv

def read_date_dirs(dpath='.',
                   reader=pd.read_csv,
                   expected_date_format=None,
                   verbose=False,
                   **kwargs):
    """Return concatenated dataframe made up of dataframes read from
    files contained in date subdirectories. 
    """
    dataframes = []
    for dname in os.listdir(dpath):
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
                    if verbose: print('Reading '+fname)
                    df = reader(os.path.join(fullpath,fname),**kwargs)
                    dataframes.append(df)
    return pd.concat(dataframes)

def plot_wind(df,
              height_name='height',
              speed_name='speed',
              direction_name='direction',
              datetime_range=(None,None)):
    """Make time-height plot of wind speed and direction"""
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
    x,y = np.meshgrid(time,height,indexing='ij')
    Nt, Nh = len(time), len(height)
    wspd = np.zeros((Nt,Nh))
    wdir = np.zeros((Nt,Nh))
    for k,h in enumerate(height):
        wspd[:,k] = dfsub.loc[dfsub['height']==h,speed_name]
        wdir[:,k] = dfsub.loc[dfsub['height']==h,direction_name]

    # make plot
    fig,ax = plt.subplots(nrows=2,sharex=True,figsize=(10,6))
    cont = ax[0].contourf(x,y,wspd, levels=np.arange(26), cmap=windspeed_colormap)
    cbar = fig.colorbar(cont, ax=ax[0], ticks=np.arange(0,26,2), label='wind speed [m/s]')
    cont = ax[1].contourf(x,y,wdir, levels=np.arange(0,361,15), cmap=winddirection_colormap)
    cbar = fig.colorbar(cont, ax=ax[1], ticks=np.arange(0,361,45), label='wind direction [deg]')

    return fig, ax
