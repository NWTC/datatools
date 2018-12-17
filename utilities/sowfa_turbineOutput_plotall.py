#!/usr/bin/env python
#
# Script to check precursor convergence and plot profiles
# written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from datatools.SOWFA.postProcessing.turbineOutput import TurbineOutput

save_dpi = 150
radial_sample_location = 0.75

turbines = TurbineOutput(ignoreHeaderNames=True,verbose=False)
if len(sys.argv) > 1:
    turbineList = [ int(iturb0) for iturb0 in sys.argv[1:] ]
else:
    turbineList = turbines.turbineList

def set_default_format(ax):
    ax[-1].set_xlabel('Time [s]')

    # set default label size
    for ax_i in ax:
        ax_i.tick_params(axis='y',labelsize='small')
    leg = ax[0].legend(fontsize='small',loc='upper left',bbox_to_anchor=(1.,1.),borderpad=0.1)

    # get rid of whitespace
    plt.tight_layout()

    # correct axes extent for outside legend
    fig.canvas.draw() # needed to get correct legend width
    legend_width = leg.get_window_extent().inverse_transformed(fig.transFigure).width
    plt.subplots_adjust(right=1.0-legend_width)


# plot rotor power, speed, torque, and pitch

rotor = turbines.readRotorOutputs()
print('Reading rotor outputs')

fig,ax = plt.subplots(nrows=4,sharex=True,figsize=(6,4))
for iturb in turbineList:
    turb = rotor[iturb]
    (turb['rotorPower']/1000).plot(ax=ax[0],label='WT'+str(iturb))
    turb['rotorSpeed'].plot(ax=ax[1])
    (turb['rotorTorque']/1000).plot(ax=ax[2])
    turb['bladePitch'].plot(ax=ax[3])
ax[0].set_ylabel('Rotor Power\n[kW]',fontsize='small')
ax[1].set_ylabel('Rotor Speed\n[rpm]',fontsize='small')
ax[2].set_ylabel('Rotor Torque\n[kN-m]',fontsize='small')
ax[3].set_ylabel('Blade Pitch\n[deg]',fontsize='small')
set_default_format(ax)
fig.savefig('turbine_outputs-rotor.png', bbox_inches='tight',dpi=save_dpi)

# plot blade 1 max velocity, angle of attack

print('Reading blade 0 outputs at 75%R')
blade0 = turbines.readBladeOutputs(bladeNum=0)
ridx = int(radial_sample_location*(turbines.Npoints-1))
rname = '{:g}'.format(radial_sample_location).replace('.','')
suffix = '_' + str(ridx)

fig,ax = plt.subplots(nrows=4,sharex=True,figsize=(6,4))
for iturb in turbineList:
    blade = blade0[iturb]
    blade['bladePointVtangential'+suffix].plot(ax=ax[0],label='WT'+str(iturb))
    blade['bladePointVaxial'+suffix].plot(ax=ax[1])
    blade['bladePointVradial'+suffix].plot(ax=ax[2])
    blade['bladePointAlpha'+suffix].plot(ax=ax[3])
ax[0].set_ylabel('Tangential\nVelocity\n[m/s]',fontsize='small')
ax[1].set_ylabel('Axial\nVelocity\n[m/s]',fontsize='small')
ax[2].set_ylabel('Radial\nVelocity\n[m/s]',fontsize='small')
ax[3].set_ylabel('Angle of Attack\n[deg]',fontsize='small')
set_default_format(ax)
fig.savefig('turbine_outputs-blade_{:s}R.png'.format(rname),
            bbox_inches='tight',dpi=save_dpi)


#---------
print('Loaded data:')
for key,val in turbines.names.items():
    print('   ',key,':',val)
plt.show()

