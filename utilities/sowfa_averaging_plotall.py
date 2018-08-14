#!/usr/bin/env python
#
# Script to check precursor convergence
# written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

from datatools.SOWFA.postProcessing.averaging import PlanarAverages

plt.style.use('seaborn-darkgrid')

heights = [90.,200, 400, 600, 800] # to sample TI history
tavg_window = 600.0 # for statistics
dt = 1.0 # for data resampling
include_SFS = True # include stresses from SGS model in TI calculation
wind_aligned = True

# reads avg.hLevelsCell with shape (NZ)
# reads avg.U_mean, ... with shape (NT, NZ)
avg = PlanarAverages( *sys.argv[1:] )

#------------------------------------------------------------------------------
#
# Histories of turbulence intensity, turbulent kinetic energy
#
avg.calculate_TI(heights, tavg_window=tavg_window, dt=dt, SFS=include_SFS)
fig,ax = avg.plot_TI_history(savefig='TIhist.png')
fig,ax = avg.plot_TKE_history(savefig='TKEhist.png')

avg.save_TI_history(prefix='TIhist')

print('last averaging time:',avg.tavg[-1])
for i,h in enumerate(heights):
    print('z= {:.1f}:'.format(h))
    print('  TIx/y/z = {:g} {:g} {:g}'.format(avg.TIx[-1,i],avg.TIy[-1,i],avg.TIz[-1,i]))
    print('  TIdir   = {:g}'.format(avg.TIdir[-1,i]))
    print('  TIxyz   = {:g}'.format(avg.TIxyz[-1,i]))
    print('  TKE     = {:g}'.format(avg.TKE[-1,i]))

# calculate rotated tensors for all times -- this is slow
#if wind_aligned:
#    avg.rotate_tensors()

#------------------------------------------------------------------------------
#
# Wind and Temperature Profiles 
#
if wind_aligned:
    fig,ax = plt.subplots(ncols=3)
    avg.plot_windspeed_profile(ax=ax[0])
    avg.plot_winddirection_profile(ax=ax[1])
    avg.plot_T_profile(ax=ax[2])
    ax[1].set_ylabel('')
    ax[2].set_ylabel('')
else:
    fig,ax = plt.subplots(ncols=2)
    avg.plot_UVW_profile(ax=ax[0])
    avg.plot_T_profile(ax=ax[1])
    ax[1].set_ylabel('')
fig.savefig('Profiles_Mean.png',bbox_inches='tight')
fig.suptitle('Resolved Mean Quantities')

#------------------------------------------------------------------------------
#
# Resolved Fluctuating Quantities
#
fig,ax = plt.subplots(ncols=2)
avg.plot_variance_profile(ax=ax[0],rotated=wind_aligned)
avg.plot_covariance_profile(ax=ax[1],rotated=wind_aligned)
ax[1].set_ylabel('')
fig.savefig('Profiles_Fluc.png',bbox_inches='tight')
fig.suptitle('Resolved Fluctuating Quantities')

#------------------------------------------------------------------------------
#
# Modeled SFS Quantities
#
fig,ax = plt.subplots(ncols=2)
avg.plot_SFS_normalstress_profile(ax=ax[0],rotated=wind_aligned)
avg.plot_SFS_shearstress_profile(ax=ax[1],rotated=wind_aligned)
ax[1].set_ylabel('')
fig.savefig('Profiles_SFS.png',bbox_inches='tight')
fig.suptitle('Sub-Filter Scale Quantities')

#------------------------------------------------------------------------------
#
# Calculated quantities
#
avg.calculate_TI_profile(tavg_window=tavg_window, dt=dt, SFS=include_SFS)

fig,ax = plt.subplots(ncols=2)
ax[0].plot(100*avg.TI_profile, avg.hLevelsCell, 'k-')
ax[1].plot(avg.TKE_profile, avg.hLevelsCell, 'k-')
ax[0].set_xlabel(r'Turbulence Intensity [%]')
ax[1].set_xlabel(r'Turbulence Kinetic Energy [m$^2$/$s^2$]')
ax[0].set_ylabel(r'Height [m]')
ax[1].set_ylabel('')
fig.savefig('Profiles_TI.png',bbox_inches='tight')

#------------------------------------------------------------------------------
#
# Save averaging data
#
avg.save_profile(fname='averagingProfiles.csv') # latest time

# pandas dataframe
avg.to_csv('averaging.csv') # all times

plt.show()

