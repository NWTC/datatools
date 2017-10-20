#!/usr/bin/env python
#
# Module for in and outputting data in the OpenFOAM timeVaryingMappedFixedValue format.
#
# Written by Eliot Quon (eliot.quon@nrel.gov) -- 2017-10/18
#
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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"""

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
{avgValue}\n\n"""

def writeBoundaryData(fname,data,
        patchName='patch',
        timeName=0,
        avgValue=None):
    """Write out a boundaryData file which should be stored in
        constant/boundarydata/patchName/timeName/fname

    NOTE: This has not been tested yet.

    Parameters
    ----------
    fname : str
        Output data file name
    data : numpy.ndarray
        Field data to be written out, with shape (3,NY,NZ) for vectors
        and shape (NY,NZ) for scalars
    patchName : str, optional
        Name of the boundary patch
    timeName : scalar or str, optional
        Name of corresponding time directory
    avgValue : scalar or list-like, optional
        To set avgValue in the data file; probably not used.

    @author: ewquon
    """
    dims = data.shape
    if len(dims) == 2:
        patchType = 'scalar'
        NY,NZ = dims
        if avgValue is None:
            avgValueStr = '0'
        else:
            avgValueStr = str(avgValue)
    elif len(dims) == 3:
        patchType = 'vector'
        NY,NZ = dims[1:]
        if avgValue is None:
            avgValueStr = '(0 0 0)'
        else:
            avgValueStr = '(' + ' '.join([str(val) for val in list(avgValue)]) + ')'
    else:
        print 'ERROR: Unexpected number of dimensions! No data written.'
        return

    with open(fname,'w') as f:
        f.write(dataheader.format(
                patchType=patchType,
                patchName=patchName,
                timeName=timeName,
                avgValue=avgValueStr))
        f.write('{:d}\n(\n'.format(NY*NZ))
        for k in range(NZ):
            for j in range(NY):
                f.write('({v[0]:f} {v[1]:f} {v[2]:f})\n'.format(v=data[:,j,k]))
        f.write(')\n')


def _getPointsFromList(ylist,zlist):
    """Detects y and z (1-D arrays) from a list of points on a
    structured grid. Makes no assumptions about the point
    ordering
    """
    ylist = np.array(ylist)
    zlist = np.array(zlist)
    N = len(zlist)
    assert(N == len(ylist))
    if zlist[1]==zlist[0]:
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
    return y,z

def readBoundaryPoints(fname,checkConst=True,tol=1e-6):
    """Returns a 2D set of points if one of the coordinates is constant
    otherwise returns a 3D set of points.
    Assumes that the points are on a structured grid.
    """
    N = None
    points = None
    iread = 0
    with open(fname,'r') as f:
        for line in f:
            if N is None:
                try:
                    N = int(line)
                    points = np.zeros((N,3))
                    print 'Reading',N,'points from',fname
                except ValueError: pass
            elif not line.strip() in ['','(',')']:
                points[iread,:] = [ float(val) for val in line.strip().strip('()').split() ]
                iread += 1
    assert(iread == N)

   #constX = np.all(points[:,0] == points[0,0])
   #constY = np.all(points[:,1] == points[0,1])
   #constZ = np.all(points[:,2] == points[0,2])
    constX = np.max(points[:,0]) - np.min(points[0,0]) < tol
    constY = np.max(points[:,1]) - np.min(points[0,1]) < tol
    constZ = np.max(points[:,2]) - np.min(points[0,2]) < tol
    print 'Constant in x/y/z :',constX,constY,constZ
    assert(constX or constY)

    if constX:
        ylist = points[:,1]
        zlist = points[:,2]
    elif constY:
        ylist = points[:,0]
        zlist = points[:,2]

    return _getPointsFromList(ylist,zlist)

def readVectorData(fname,NY=None,NZ=None,order='F'):
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
                    print 'Reading',N,'vectors from',fname
                except ValueError: pass
            elif not line.strip() in ['','(',')']:
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


def readScalarData(fname,NY=None,NZ=None,order='F'):
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
                        print 'Reading',N,'scalars from',fname
                except ValueError: pass
            elif not line.strip() in ['','(',')']:
                data[iread] = float(line)
                iread += 1
    assert(iread == N)

    if (NY is not None) and (NZ is not None):
        scalarField = data.reshape((NY,NZ),order=order)
    else:
        scalarField = data

    return scalarField

