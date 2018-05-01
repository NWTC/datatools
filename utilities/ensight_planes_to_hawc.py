#!/usr/bin/env python
#
# Script to convert "array" sampling planes from SOWFA (written in
# Ensight format) and into a FAST turbulence box (written in binary
# HAWC full-field file format)
#
# Written by Eliot Quon (eliot.quon@nrel.gov)
#
# SAMPLE USAGE:
# To process ./*/inflowPlane_03km_U.000.U:
#   ensight_planes_to_hawc.py 'inflowPlane_03km_U'
#
# Note: Alternative values for reference velocity and height (uref
# and zref) may be specified but should not be modify the velocity
# field. They are subtracted out by the this script only to be added
# back on by FAST (based on what is specified in the InflowWind input
# file).
#
from __future__ import print_function
import numpy as np

from datatools.dataloaders import foam_ensight_array
from datatools.FAST.InflowWind import input_template
from datatools.binario import binaryfile


def generate_inflow(prefix,
        uref=8.0,zref=90.0,
        ufile='u.bin',vfile='v.bin',wfile='w.bin',
        inflowfile='InflowWind_from_SOWFA.dat'):
    """Writes out one binary file for each wind component in the HAWC
    format as described in the InflowWind manual, in addition to an
    InflowWind input file"""

    inflow = foam_ensight_array('.', prefix=prefix,
                                npzdata=prefix+'.npz') # auto-detect NX,NY,NZ

    # time series detected from directory names
    t = np.array(inflow.ts.outputTimes)

    # assume flow is in x direction
    X,Y,Z,U = inflow.sliceI() # return arrays with shape (NY,NZ)
                              # or in the case of U: (Ntimes,NY,NZ,3)
    assert(np.min(X) == np.max(X)) # plane is at constant x

    # calculate turbulence box description
    nt = inflow.ts.Ntimes
    ny = inflow.NY
    nz = inflow.NZ
    y = Y[:,0]
    z = Z[0,:]
    print('x :',nt,(t-t[0])*uref)
    print('y :',ny,y)
    print('z :',nz,z)
    dx = uref*(t[1]-t[0])
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    U[:,:,:,0] -= uref # InflowWind will add this back to the x-component
    with binaryfile(ufile,'w') as f:
        # last plane of turbulence box enters rotor first, and corresponds to
        # the first time snapshot
        for i in range(nt): # indexing goes nx, nx-1, ... 1
            for j in range(ny)[::-1]: # indexing goes ny, ny-1, ... 1
                f.write_float(U[i,j,:,0]) # indexing goes 1, 2, ... nz
    print('Wrote binary',ufile)

    with binaryfile(vfile,'w') as f:
        # last plane of turbulence box enters rotor first, and corresponds to
        # the first time snapshot
        for i in range(nt): # indexing goes nx, nx-1, ... 1
            for j in range(ny)[::-1]: # indexing goes ny, ny-1, ... 1
                f.write_float(U[i,j,:,1]) # indexing goes 1, 2, ... nz
    print('Wrote binary',vfile)

    with binaryfile(wfile,'w') as f:
        # last plane of turbulence box enters rotor first, and corresponds to
        # the first time snapshot
        for i in range(nt): # indexing goes nx, nx-1, ... 1
            for j in range(ny)[::-1]: # indexing goes ny, ny-1, ... 1
                f.write_float(U[i,j,:,2]) # indexing goes 1, 2, ... nz
    print('Wrote binary',wfile)

    with open(inflowfile,'w') as f:
        f.write(
            input_template.format(
                WindType=5,
                RefHt=zref,
                URef=uref,
                hawc_ufile=ufile,
                hawc_vfile=vfile,
                hawc_wfile=wfile,
                nx=nt, ny=ny, nz=nz,
                dx=dx, dy=dy, dz=dz
        ))
    print('Wrote',inflowfile)


#===============================================================================
if __name__ == '__main__':
    import sys
    prefix = sys.argv[1]
    generate_inflow(prefix)

