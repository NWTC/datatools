#!/usr/bin/env python
#
# Module for in and outputting data in the OpenFOAM timeVaryingMappedFixedValue format.
#
# Written by Eliot Quon (eliot.quon@nrel.gov) -- 2017-10/18
#
from __future__ import print_function
import numpy as np

pointsheader = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.x                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       vectorField;
    location    "constant/boundaryData/{patchName}";
    object      points;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{N}
("""

dataheader = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.x                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {patchType}AverageField;
    location    "constant/boundaryData/{patchName}/{timeName}";
    object      values;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Average
{avgValue}

{N}
("""

def write_points(fname,x,y,z,patchName='patch'):
    """Write out a points file which should be stored in
        constant/boundaryData/patchName/points
    """
    N = len(x)
    assert(N == len(y) == len(z))
#    with open(fname,'w') as f:
#        f.write(pointsheader.format(patchName=patchName,N=N))
#        f.write('{:d}\n(\n'.format(N))
#        for i in range(N):
#            f.write('({:f} {:f} {:f})\n'.format(x[i],y[i],z[i]))
#        f.write(')\n')
    np.savetxt(fname,
               np.stack((x,y,z)).T, fmt='(%f %f %f)',
               header=pointsheader.format(patchName=patchName,N=N),
               footer=')',
               comments='')

def write_data(fname,
               data,
               patchName='patch',
               timeName=0,
               avgValue=None):
    """Write out a boundaryData file which should be stored in
        constant/boundarydata/patchName/timeName/fieldName

    Parameters
    ----------
    fname : str
        Output data file name
    data : numpy.ndarray
        Field data to be written out, with shape (3,N) for vectors
        and shape (N) for scalars; 2-D or 3-D data should be flattened
        beforehand.
    patchName : str, optional
        Name of the boundary patch
    timeName : scalar or str, optional
        Name of corresponding time directory
    avgValue : scalar or list-like, optional
        To set avgValue in the data file; probably not used.

    @author: ewquon
    """
    dims = data.shape
    N = dims[-1]
    if len(dims) == 1:
        patchType = 'scalar'
        if avgValue is None:
            avgValueStr = '0'
        else:
            avgValueStr = str(avgValue)
    elif len(dims) == 2:
        patchType = 'vector'
        assert(dims[0] == 3)
        if avgValue is None:
            avgValueStr = '(0 0 0)'
        else:
            avgValueStr = '(' + ' '.join([str(val) for val in list(avgValue)]) + ')'
    else:
        print('ERROR: Unexpected number of dimensions! No data written.')
        return

#    with open(fname,'w') as f:
#        f.write(dataheader.format(patchType=patchType,
#                                  patchName=patchName,
#                                  timeName=timeName,
#                                  avgValue=avgValueStr,
#                                  N=N))
#        f.write('{:d}\n(\n'.format(N))
#        if patchType == 'vector':
#            for i in range(N):
#                f.write('({v[0]:g} {v[1]:g} {v[2]:g})\n'.format(v=data[:,i]))
#        elif patchType == 'scalar':
#            for i in range(N):
#                f.write('{v:g}\n'.format(v=data[i]))
#        f.write(')\n')

    headerstr = dataheader.format(patchType=patchType,
                                  patchName=patchName,
                                  timeName=timeName,
                                  avgValue=avgValueStr,
                                  N=N)
    if patchType == 'vector':
        np.savetxt(fname,
                   data.T, fmt='(%g %g %g)',
                   header=headerstr, footer=')',
                   comments='')
    elif patchType == 'scalar':
        np.savetxt(fname,
                   data.reshape((N,1)), fmt='%g',
                   header=headerstr, footer=')',
                   comments='')

def _get_unique_points_from_list(ylist,zlist,NY=None,NZ=None,order='F'):
    """Detects y and z (1-D arrays) from a list of points on a
    structured grid. Makes no assumptions about the point
    ordering
    """
    ylist = np.array(ylist)
    zlist = np.array(zlist)
    N = len(zlist)
    assert(N == len(ylist))
    if (NY is not None) and (NZ is not None):
        # use specified plane dimensions
        assert(NY*NZ == N)
        y = ylist.reshape((NY,NZ))[:,0]
    elif zlist[1]==zlist[0]:
        # y changes faster, F-ordering
        NY = np.nonzero(zlist > zlist[0])[0][0]
        NZ = int(N / NY)
        assert(NY*NZ == N)
        y = ylist[:NY]
        z = zlist.reshape((NY,NZ),order='F')[0,:]
    elif ylist[1]==ylist[0]:
        # z changes faster, C-ordering
        NZ = np.nonzero(ylist > ylist[0])[0][0]
        NY = int(N / NZ)
        assert(NY*NZ == N)
        z = zlist[:NZ]
        y = ylist.reshape((NY,NZ),order='C')[:,0]
    else:
        print('Unrecognized point distribution')
        print('"y" :',len(ylist),ylist)
        print('"z" :',len(zlist),zlist)
        return ylist,zlist,False
    return y,z,True

def read_boundary_points(fname,checkConst=True,tol=1e-6,**kwargs):
    """Returns a 2D set of points if one of the coordinates is constant
    otherwise returns a 3D set of points.
    Assumes that the points are on a structured grid.
    """
    N = None
    points = None
    iread = 0
    with open(fname,'r') as f:
        while N is None:
            try:
                N = int(f.readline())
            except ValueError: pass
            else:
                points = np.zeros((N,3))
                print('Reading',N,'points from',fname)
        for line in f:
            line = line[:line.find('\\')].strip()
            try:
                points[iread,:] = [ float(val) for val in line[1:-1].split() ]
            except (ValueError, IndexError): pass
            else:
                iread += 1
    assert(iread == N)

   #constX = np.all(points[:,0] == points[0,0])
   #constY = np.all(points[:,1] == points[0,1])
   #constZ = np.all(points[:,2] == points[0,2])
    constX = np.max(points[:,0]) - np.min(points[0,0]) < tol
    constY = np.max(points[:,1]) - np.min(points[0,1]) < tol
    constZ = np.max(points[:,2]) - np.min(points[0,2]) < tol
    print('Constant in x/y/z :',constX,constY,constZ)
    if not (constX or constY):
        print('Warning: boundary is not constant in X or Y?')

    if constX:
        ylist = points[:,1]
        zlist = points[:,2]
    elif constY:
        ylist = points[:,0]
        zlist = points[:,2]
    elif constZ:
        ylist = points[:,0]
        zlist = points[:,1]
    else:
        print('Unexpected boundary orientation, returning full list of points')
        return points

    return _get_unique_points_from_list(ylist,zlist,**kwargs)

def read_vector_data(fname,NY=None,NZ=None,order='F'):
    N = None
    data = None
    iread = 0
    with open(fname,'r') as f:
        for line in f:
            if N is None:
                try:
                    N = int(line)
                    if (NY is not None) and (NZ is not None):
                        if not N == NY*NZ:
                            NY = None
                            NZ = None
                    data = np.zeros((N,3))
                    print('Reading',N,'vectors from',fname)
                except ValueError: pass
            elif not line.strip() in ['','(',')',';'] \
                    and not line.strip().startswith('//'):
                data[iread,:] = [ float(val) for val in line.strip().strip('()').split() ]
                iread += 1
    assert(iread == N)

    if (NY is not None) and (NZ is not None):
        vectorField = np.zeros((3,NY,NZ))
        for i in range(3):
            vectorField[i,:,:] = data[:,i].reshape((NY,NZ),order=order)
    else:
        vectorField = data.T

    return vectorField


def read_scalar_data(fname,NY=None,NZ=None,order='F'):
    N = None
    data = None
    iread = 0
    with open(fname,'r') as f:
        for line in f:
            if (N is None) or N < 0:
                try:
                    if N is None: 
                        avgval = float(line)
                        N = -1 # skip first scalar, which is the average field value (not used)
                    else:
                        assert(N < 0)
                        N = int(line) # now read the number of points
                        if (NY is not None) and (NZ is not None):
                            if not N == NY*NZ:
                                NY = None
                                NZ = None
                        data = np.zeros(N)
                        print('Reading',N,'scalars from',fname)
                except ValueError: pass
            elif not line.strip() in ['','(',')',';'] \
                    and not line.strip().startswith('//'):
                data[iread] = float(line)
                iread += 1
    assert(iread == N)

    if (NY is not None) and (NZ is not None):
        scalarField = data.reshape((NY,NZ),order=order)
    else:
        scalarField = data

    return scalarField

