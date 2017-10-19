#!/usr/bin/env python
#
# Module for in and outputting data in the OpenFOAM timeVaryingMappedFixedValue format.
#
# Written by Eliot Quon (eliot.quon@nrel.gov) -- 2017-10/18
#

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

