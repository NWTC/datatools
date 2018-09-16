#
# Module to handle SOWFA boundary data that belong in
#   casedir/constant/boundaryData
#
# Written by Eliot Quon (eliot.quon@nrel.gov)
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
    if len(x.shape) > 1:
        x = x.ravel()
        y = y.ravel()
        z = z.ravel()
        N = len(x)
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
