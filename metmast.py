"""
Data readers for meteorological masts
=====================================
Written by Eliot Quon (eliot.quon@nrel.gov)

This is a collection of readers to be used with the NWTC datatools.wfip2
module for processing WFIP2 data downloaded from the A2e Data Archive
and Portal (DAP). No effort is made to standardize the names of the 
dataframe columns and the original data headers are retained wherever
possible.

"""
import os
import numpy as np
import pandas as pd


def load_sonic_data(fname,height=None,calculate_speed_direction=True):
    df = pd.read_csv(fname, parse_dates=['datetime'])
    df.set_index('datetime', inplace=True)
    if height is not None:
        df['z'] = height
    df['speed'] = np.sqrt(df['u']**2 + df['v']**2 + df['w']**2)
    df['direction'] = 180.0/np.pi*np.arctan2(-df['u'],-df['v'])
    df.loc[df['direction'] < 0, 'direction'] += 360.0
    return df


def upsample_with_index(df_15min,index,columns,endtime=None):
    """Given a dataframe sampled at a lower frequency (e.g., 15 min), 
    create a new dataframe with the specified index--presumably at a 
    higher sampling frequency (e.g., 1 s).

    The speciifed data columns are copied from the low-frequency data-
    frame to the high-frequency dataframe.

    Intermediate values are interpolated.
    """
    # create dataframe with indices from 2016-11-01 00:00:00 to
    # 2016-11-30 23:59:59 (1-s intervals)
    df = pd.DataFrame(index=index)
    # set values at 15-min intervals
    df[columns] = df_15min[columns]
    # fill nans
    if (endtime is not None) and (df_15min.index[-1] == endtime):
        # add row at 2016-12-01 00:00:00 since if available
        df = df.append(df_15min.iloc[-1])
        df = df.interpolate(method='slinear')
        df = df.iloc[:-1]
    else:
        df = df.interpolate(method='slinear')
    assert(np.all(df.index == index))
    return df


cols_no_time = ['u','v','w','Tv','qc']

def firstrow_freq_reader(fpath,
                         prefix='met.z09.b0.',suffix='.csv',
                         columns=cols_no_time,
                         freq=20.0,
                         na_values=99.99,
                         qcflag=0):
    """Given a sampling frequency, freq [Hz], take every int(freq) rows
    to effectively sample at 1 Hz. Return only data that has been flagged
    with qcflag (if qcflag is not None).
    """
    fname = os.path.split(fpath)[-1]
    assert(fname.startswith(prefix) and fname.endswith(suffix))
    datetime_str = fname[len(prefix):-len(suffix)]
    datetime0 = pd.to_datetime(datetime_str, format='%Y%m%d.%H%M%S')
    df = pd.read_csv(fpath, names=columns, na_values=na_values)
    N = int(freq)
    df = df.iloc[::N,:]
    df['datetime'] = datetime0 + pd.timedelta_range(start='0s',periods=len(df),freq='1s')
    df['Tv'] += 273.15
    df = df.set_index('datetime')
    if qcflag is None:
        return df
    else:
        columns.remove('qc') 
        return df.loc[df['qc']==qcflag, columns]


cols_ANL = ['X','Y','Z','T','time']

def firstrow_reader(fpath,prefix='met.z28.b0.',suffix='.sonic80ms.txt',checkdate=True):
    """Given times in HH:MM:SS format with duplicate timestamps,
    save the first row to effectively resample at 1 Hz
    """
    fname = os.path.split(fpath)[-1]
    assert(fname.startswith(prefix) and fname.endswith(suffix))
    datetime_str = fname[len(prefix):-len(suffix)]
    datetime0 = pd.to_datetime(datetime_str, format='%Y%m%d.%H%M%S')
    lasttime = ''
    data = []
    with open(fpath,'r') as f:
        for line in f:
            row = line.split()
            if not row[-1] == lasttime:
                data.append(row)
                lasttime = row[-1]
    df = pd.DataFrame(data=data,columns=cols_ANL)
    df['datetime'] = datetime0.date() + pd.to_timedelta(df['time'])
    if checkdate:
        assert(df.iloc[0]['datetime'] == datetime0)
    elif not df.iloc[0]['datetime'] == datetime0:
        print('Missing data? Range: {} {} ({})'.format(df.iloc[0]['datetime'],
                                                       df.iloc[-1]['datetime'],
                                                       fname))
    df['u'] = pd.to_numeric(df['Y']) / 100.0  # [m/s]
    df['v'] = pd.to_numeric(df['X']) / 100.0  # [m/s]
    df['w'] = pd.to_numeric(df['Z']) / 100.0  # [m/s]
    df['Tv'] = pd.to_numeric(df['T']) / 100.0 + 273.15  # [K]
    df = df.set_index('datetime')
    return df[['u','v','w','Tv']]


# def RMYoung_height_wrapper(dpath,prefix,heights,*args):
#     """call read_RMYoung_sonic (raw data) for different heights"""
#     dflist = []
#     for h in heights:
#         df = wfip2.read_date_dirs(dpath=dpath,
#                                   reader=read_RMYoung81000_sonic,
#                                   ext='.son{:02d}m.dat'.format(int(h)))
#         df['z'] = h
#         dflist.append(df)
#     return pd.concat(dflist)

cols_RMY81000 = ['Y','m','d','H','M','S','ms','record','sonic_u','sonic_v','w','sonic_temperature']

def RMYoung81000_sonic(*args):
    """E.g., to read UND RM Young 81000 sonics"""
    dflist = []
    for fname in args:
        df = pd.read_csv(fname,names=cols_RMY81000)
        #df['us'] = df['ms'] * 1000.0
        df = df.loc[df['ms']==0]
        df['datetime'] = pd.to_datetime(
                df.Y*10000000000 + df.m*100000000 + df.d*1000000 +
                df.H*10000 + df.M*100 + df.S,
                format='%Y%m%d%H%M%S')
        #df['datetime'] = str(df['Y'])+str(df['m'])+str(df['d']) + str(df['H'])+str(df['M'])+str(df['S'])
        #df['datetime'] = pd.to_datetime(df['datetime'])
        df['u'] = -df['sonic_u']
        df['v'] = -df['sonic_v']
        df['th'] = df['sonic_temperature'] + 273.15
        dflist.append( df[['datetime','u','v','w','th']] )
    return pd.concat(dflist).set_index('datetime')


cols_RMY05106 = ['id',
             'year','julian_day','HHMM',
             'S', # mean horizontal wind speed [m/s]
             'U', # resultant mean wind speed [m/s]
             'dir', # resultant mean of wind direction == arctan(Ueast/Unorth) [deg]
             'dir_std', # stdev of wind direction [deg]
             'T', # air temperature [C]
             'RH', # relative humidity [%]
             'P', # barometric pressure (not corrected for sea level) [mb]
             'E', # downwelling SW solar radiation (400-1100nm) [W/m^2]
             'T10X', # datalogger temperature [C]
             'p10X', # datalogger power [V]
            ]

def RMYoung05106_sonic(*args):
    """E.g., to read RM Young 05106 sonics"""
    dflist = []
    for fname in args:
        df = pd.read_csv(fname,names=cols_RMY05106)
        df['datetime'] = pd.to_datetime(
                df.year*10000000 + df.julian_day*10000 + df.HHMM,
                format='%Y%j%H%M')
        ang = np.radians(270 - df['dir'])
        df['u'] = df['U'] * np.cos(ang)
        df['v'] = df['U'] * np.sin(ang)
        df['th'] = df['T'] + 273.15
        df = df.rename({'P':'p'},axis=1)
        dflist.append( df[['datetime','u','v','th','p','E']] )
    return pd.concat(dflist).set_index('datetime')


# def eddypro_height_wrapper(dpath,prefix,heights,*args,**kwargs):
#     """call read_eddypro, i.e., RM Young sonic data postprocessed using
#     EddyPro, for different heights
#     """
#     dflist = []
#     for h in heights:
#         df = wfip2.read_dir(dpath=dpath,
#                             reader=read_eddypro,
#                             prefix=prefix,
#                             ext='.son{:02d}m.full_output.csv'.format(int(h)),
#                             **kwargs)
#         df['z'] = h
#         dflist.append(df)
#     return pd.concat(dflist)

def EddyPro_full_output(*args,**kwargs):
    """read EddyPro "full output" csv files"""
    columns = kwargs.get('columns',None)
    rename = kwargs.get('rename',False)
    verbose = kwargs.get('verbose',False)
    if (columns is not None) and ('datetime' not in columns):
        columns = ['datetime'] + columns
    dflist = []
    for fname in args:
        if verbose:
            print('reading',fname)
        with open(fname,'r') as f:
            f.readline()
            header = f.readline().strip().split(',')
            f.readline() # units
            df = pd.read_csv(f,names=header,parse_dates={'datetime':['date','time']},na_values=-9999)
            if columns is not None:
                df = df[columns]
            if rename is not False:
                if rename is True:
                    renamedict = {
                        'wind_speed':'windspeed',
                        'wind_dir':'winddir',
                        'sonic_temperature':'temperature',
                    }
                else:
                    renamedict = rename
                df = df.rename(columns=renamedict)
            dflist.append(df)
    return pd.concat(dflist).set_index('datetime')


