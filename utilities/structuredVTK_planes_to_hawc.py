#!/usr/bin/env python
#
# Script to convert structured VTK data from SOWFA into a FAST
# turbulence box (written in binary HAWC full-field file format)
#
# Written by Eliot Quon (eliot.quon@nrel.gov)
# Modified by Matt Churchfield (matt.churchfield@nrel.gov)
#
from __future__ import print_function
import numpy as np

from datatools.dataloaders import foam_structuredVTK_array
from datatools.FAST.InflowWind import input_template
from datatools.binario import binaryfile


def generate_inflow(datadir,prefix,
        uref=8.0,zref=90.0,
        tstart=None,tend=None,
        ufile='u.bin',vfile='v.bin',wfile='w.bin',
        inflowfile='InflowWind_from_SOWFA.dat'):
    """Writes out one binary file for each wind component in the HAWC
    format as described in the InflowWind manual, in addition to an
    InflowWind input file"""

    inflow = foam_structuredVTK_array(datadir, prefix=prefix,
                                      tstart=tstart, tend=tend,
                                      npzdata=prefix+'.npz')

    # time series detected from directory names
    t = np.array(inflow.ts.times)

    # assume flow is in x direction
    X,Y,Z,U = inflow.sliceI() # return arrays with shape (NY,NZ)
                              # or in the case of U: (Ntimes,NY,NZ,3)
    assert(np.min(X) == np.max(X)) # plane is at constant x
    
    #Create a slice of the time record
#    indices = np.nonzero((t >= tstart) & (t <= tend))[0]
#    t = t[indices]
#    print 'selected time range:',np.min(t),np.max(t)
#    U = U[indices,:,:,:]

    # calculate turbulence box description
#    nx = len(t) # inflow.ts.Ntimes
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



    jProbe = 160
    kProbe = [0, 10, 30, 70, 100]
    pf = open('probeFileU.dat','w')
    for i in range(nt):
       pf.write(str(t[i]) + ' ')
       for k in range(len(kProbe)):
           if (k < len(kProbe)-1):
               pf.write(str(U[i,jProbe,kProbe[k],0]) + ' ')
           else:
               pf.write(str(U[i,jProbe,kProbe[k],0]) + '\n')
    pf.close()
    print('Wrote u probe file')

    pf = open('probeFileV.dat','w')
    for i in range(nt):
       pf.write(str(t[i]) + ' ')
       for k in range(len(kProbe)):
           if (k < len(kProbe)-1):
               pf.write(str(U[i,jProbe,kProbe[k],1]) + ' ')
           else:
               pf.write(str(U[i,jProbe,kProbe[k],1]) + '\n')
    pf.close()
    print('Wrote v probe file')

    pf = open('probeFileW.dat','w')
    for i in range(nt):
       pf.write(str(t[i]) + ' ')
       for k in range(len(kProbe)):
           if (k < len(kProbe)-1):
               pf.write(str(U[i,jProbe,kProbe[k],2]) + ' ')
           else:
               pf.write(str(U[i,jProbe,kProbe[k],2]) + '\n')
    pf.close()
    print('Wrote w probe file')




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

