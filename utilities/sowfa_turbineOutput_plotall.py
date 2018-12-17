#!/usr/bin/env python
#
# Script to check precursor convergence and plot profiles
# written by Eliot Quon (eliot.quon@nrel.gov)
#
import sys
import numpy as np
import matplotlib.pyplot as plt
from datatools.SOWFA.postProcessing.turbineOutput import TurbineOutput

td = TurbineOutput(ignoreHeaderNames=True)
if len(sys.argv) > 1:
    turbineList = [ int(iturb0) for iturb0 in sys.argv[1:] ]
else:
    turbineList = td.turbineList

outputs = td.readRotorOutputs()

fig,ax = plt.subplots(nrows=4,sharex=True,figsize=(6,4))
for iturb in turbineList:
    turb = outputs[iturb]
    (turb['rotorPower']/1000).plot(ax=ax[0],label='WT'+str(iturb))
    turb['rotorSpeed'].plot(ax=ax[1])
    (turb['rotorTorque']/1000).plot(ax=ax[2])
    turb['bladePitch'].plot(ax=ax[3])
ax[0].set_ylabel('Rotor\nPower\n[kW]')
ax[1].set_ylabel('Rotor\nSpeed\n[rpm]')
ax[2].set_ylabel('Rotor\nTorque\n[kN-m]')
ax[3].set_ylabel('Blade\nPitch\n[deg]')
ax[-1].set_xlabel('Time [s]')
leg = ax[0].legend(fontsize='small',loc='upper left',bbox_to_anchor=(1,1))

plt.tight_layout()
fig.savefig('turbine_outputs.png',bbox_extra_artists=(leg,),bbox_inches='tight')

plt.show()
