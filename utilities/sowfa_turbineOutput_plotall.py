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

fig,ax = plt.subplots(nrows=3,sharex=True,figsize=(6,4))
for iturb in turbineList:
    turb = outputs[iturb]
    P = turb['rotorPower'] / 1000.
    P.plot(ax=ax[0],label='Turbine'+str(iturb))
    turb['rotorSpeed'].plot(ax=ax[1])
    turb['bladePitch'].plot(ax=ax[2])

ax[0].set_ylabel('Rotor Power [kW]')
ax[1].set_ylabel('Rotor Speed [rpm]')
ax[2].set_ylabel('Blade Pitch [deg]')
ax[1].set_xlabel('Time [s]')
leg = ax[0].legend(fontsize='small',loc='upper left',bbox_to_anchor=(1,1))

plt.tight_layout()
fig.savefig('turbine_outputs.png',bbox_extra_artists=(leg,),bbox_inches='tight')

plt.show()
