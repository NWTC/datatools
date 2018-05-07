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

heights = [90.,200, 400, 600, 800]

# read avg.hLevelsCell with shape (NZ)
# read avg.U_mean, ... with shape (NT, NZ)
avg = PlanarAverages( *sys.argv[1:] )

#------------------------------------------------------------------------------
avg.calculate_TI(heights, tavg_window=600.0, dt=1.0, SFS=True)
fig,ax = avg.plot_TI_history(savefig='TIhist.png')

avg.save_TI_history(prefix='TIhist')

print('last averaging time:',avg.tavg[-1])
for i,h in enumerate(heights):
    print('z=',h,':')
    print('  TIx/y/z =',avg.TIx[-1,i],avg.TIy[-1,i],avg.TIz[-1,i])
    print('  TIdir   =',avg.TIdir[-1,i])
    print('  TIxyz   =',avg.TIxyz[-1,i],' ( TKE=',avg.TKE[-1,i],')')

#------------------------------------------------------------------------------
#
# Velocity and Temperature Profiles 
#
fig,ax = plt.subplots(ncols=2)
avg.plot_UVW_profile(ax=ax[0])
avg.plot_T_profile(ax=ax[1])
ax[1].set_ylabel('')
fig.suptitle('Resolved Mean Quantities')
fig.savefig('Profiles_Mean.png')

#------------------------------------------------------------------------------
#
# Wind speed and direction
#
fig,ax = plt.subplots(ncols=2)
avg.plot_windspeed_profile(ax=ax[0])
avg.plot_winddirection_profile(ax=ax[1])
ax[1].set_ylabel('')
fig.suptitle('Resolved Mean Wind Profiles')
fig.savefig('Profiles_MeanUdir.png')

#------------------------------------------------------------------------------
#
# Resolved Fluctuating Quantities
#
fig,ax = plt.subplots(ncols=2)
avg.plot_variance_profile(ax=ax[0])
avg.plot_covariance_profile(ax=ax[1])
ax[1].set_ylabel('')
fig.suptitle('Resolved Fluctuating Quantities')
fig.savefig('Profiles_Fluc.png')

#------------------------------------------------------------------------------
#
# Modeled SFS Quantities
#
fig,ax = plt.subplots(ncols=2)
avg.plot_SFS_normalstress_profile(ax=ax[0])
avg.plot_SFS_shearstress_profile(ax=ax[1])
ax[1].set_ylabel('')
fig.suptitle('Sub-Filter Scale Quantities')
fig.savefig('Profiles_SFS.png')

#------------------------------------------------------------------------------

avg.save_profile(fname='averagingProfiles.csv')

plt.show()

