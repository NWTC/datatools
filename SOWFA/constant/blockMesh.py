#
# Module for assisting with blockMesh grid generation
#
# Written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import numpy as np

def grow_mesh(N, d0, r):
    """Calculate the cell spacings given a number of cells (N), an 
    initial spacing (d0), and a normal growth rate (r).

    Spacings are calculated as follows:
        d_i = d_0 * r^i,  for i in 0 .. N-1
    """
    spacings = d0 * np.array([r**i for i in range(N)])
    print('spacings : {:g} .. {:g}'.format(d0,spacings[-1]))
    return spacings

def points_from(spacings, x0=0.0):
    x = x0 * np.ones((len(spacings)+1,))
    x[1:] += np.cumsum(spacings)
    return x

def spacings(L, N, ratio=1):
    """Calculate the blockmesh spacings that would result given a length
    (L), number of cells (N), and simpleGrading ratio.
    """
    r = np.exp(np.log(ratio) / (N-1))
    print('growth rate = ',r)
    d0 = L / np.sum([r**i for i in range(N)])
    return grow_mesh(N, d0, r)

def start_end(L, N, ratio):
    r = np.exp(np.log(ratio) / (N-1))
    d0 = L / np.sum([r**i for i in range(N)])
    return d0, ratio*d0

def estimate(d, d0, r, L=None, round_to=10., output='blockMesh'):
    """Estimate the blockmesh inputs that would give the desired final
    point spacing (d) over a specified distance (L) with approximately
    the normal growth rate (r). The actual growth rate will not exactly
    match the input unless d_n == d0 * r^(N-1). 
    
    The initial spacing (d0) is used to estimate the number of points
    within the region, and is then adjusted to give the desired final
    spacing (d) over distance (L).

    If L is not specified, then a reasonable estimate of blockMesh
    inputs is provided with L rounded to a "nice" number.

    Returns estimated N, d0, r
    """
    # first estimate the number of cells to get closed to the initial
    # and final spacings given a particular growth rate
    N_est = np.log(d/d0) / np.log(r) + 1
    N = int(np.ceil(N_est))
    print('estimated number of cells is {:g} ~= {:d}'.format(N_est, N))
    L1 = np.sum(grow_mesh(N,d0,r))
    print('resultant layer height : ',L1)
    # we assume that the estimated number of cells is pretty close to
    # what it should be and leave this fixed; since the growth rate
    # r = f(d0,d1,N), we should update it assuming the other parameters
    # are also fixed
    r_approx = np.exp(np.log(d/d0) / (N-1))
    print('actual growth rate: ',r_approx)
    # we're basically done at this point, but the resulting distance 
    # spanned by the points probably isn't a nice floating point number
    if L is None:
        L = np.round(L1/round_to) * round_to
    print('calculating d0 for L =',L)
    d0_approx = L / np.sum([r_approx**i for i in range(N)])
    print('adjusted initial spacing :',d0_approx)
    if output=='growth':
        return N, d0_approx, r_approx
    else:
        return L, N, r_approx**(N-1)


blockMeshDict_header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  1.6                                   |
|   \\\\  /    A nd           | Web:      http://www.OpenFOAM.org               |
|    \\\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include        "../../setUp"

convertToMeters 1.0;

"""

patch_def = """    {name:s}
    {{
        type patch;
        faces
        (
{faceslist:s}
        );
    }}
"""

blockMeshDict_footer = """
// ************************************************************************* //"""


class BlockMeshDict(object):
    """Used to manage mesh regions within a blockMesh definition. At the
    moment, only vertically layered regions are supported.
    """

    def __init__(self,xMin,xMax,yMin,yMax,zMin,zMax,dx,dy,dz): 
        """Initialize domain with specified overall dimensions
        (Lx,Ly,Lz) with finest-level grid spacings (dx,dy,dz). Actual
        domain extents will be greater than or equal to the inputs,
        depending on whether the input spacings exactly fill up the
        domain.
        """
        self.xMin = xMin
        self.yMin = yMin
        self.zMin = zMin
        self.xMax = xMax
        self.yMax = yMax
        self.zMax = zMax
        self.dx0 = dx
        self.dy0 = dy
        self.dz0 = dz

        self.vertex0 = 0
        self.blocks = []

        # descriptors of each region
        self.Nx = []
        self.Ny = []
        self.Nz = []
        self.z0 = []
        self.z1 = []
        self.simpleGradingZ = []

    def __repr__(self):
        Nlayers = len(self.Nx)
        s = 'Mesh bounding box corners: ({:g} {:g} {:g}) ({:g} {:g} {:g})\n'.format(
                self.xMin, self.yMin, self.zMin,
                self.xMax, self.yMax, self.zMax)
        if Nlayers == 0:
            s += '  no layers added; call generate_uniform_grid or generate_region'
        for i in range(Nlayers):
            dx = (self.xMax - self.xMin) / self.Nx[i]
            dy = (self.yMax - self.yMin) / self.Ny[i]
            dz = (self.z1[i] - self.z0[i]) / self.Nz[i]
            s+= '  layer {:d} : z=({:g} {:g}) N=({:d} {:d} {:d}) spacings=({:g} {:g} {:g})\n'.format(
                    i, self.z0[i], self.z1[i],
                    self.Nx[i], self.Ny[i], self.Nz[i],
                    dx, dy, dz)
        return s
            

    def generate_uniform_grid(self):
        self.generate_region(self.zMin, self.zMax, self.dx0, self.dy0, self.dz0,
                             add=True)

    def generate_region(self,z0,z1,dx,dy,dz0,dz1=None,r=1,add=False):
        """Generate new region between height z0 and z1, with grid
        spacings dx, dy, and dz.

        Optionally, a growth rate r > 1 may be specified, in which case
        the target final grid spacing (dz1) should be specified.

        The new region is added to the domain until the add=True is
        specified.
        """
        Nx = int(np.ceil((self.xMax-self.xMin)/dx))
        x0 = self.xMin
        x1 = x0 + Nx*dx
        if x1 > self.xMax:
            self.xMax = x1
            print('Updated xMax to',self.xMax)
        Ny = int(np.ceil((self.yMax-self.yMin)/dy))
        y0 = self.yMin
        y1 = y0 + Ny*dy
        if y1 > self.yMax:
            self.yMax = y1
            print('Updated yMax to',self.yMax)

        if r==1:
            Nz = int(np.ceil((z1-z0)/dz0))
            ratio = 1.0
        else:
            assert(dz1 is not None)
            Lz,Nz,ratio = estimate(dz1,dz0,r,L=(z1-z0))
            assert(Lz == z1-z0)
        
        if add:
            self.Nx.append(Nx)
            self.Ny.append(Ny)
            self.Nz.append(Nz)
            self.z0.append(z0)
            self.z1.append(z1)
            self.simpleGradingZ.append(ratio)
            hexblock = hex(self.vertex0, N=(Nx,Ny,Nz),
                           simpleGrading=(1,1,ratio))
            self.blocks.append(hexblock)
            self.vertex0 += 8
            print('Added region',len(self.Nx))

    def _vertices(self):
        s = 'vertices\n('
        for i in range(self.Nlayers):
            s += """
    ( $xMin   $yMin   $zMin{i:d} )
    ( $xMax   $yMin   $zMin{i:d} )
    ( $xMax   $yMax   $zMin{i:d} )
    ( $xMin   $yMax   $zMin{i:d} )
    ( $xMin   $yMin   $zMax{i:d} )
    ( $xMax   $yMin   $zMax{i:d} )
    ( $xMax   $yMax   $zMax{i:d} )
    ( $xMin   $yMax   $zMax{i:d} )""".format(i=i)
        s += '\n);\n\n'
        return s

    def _blocks(self):
        s = 'blocks\n(\n'
        for block in self.blocks:
            s += '    ' + block.write() + '\n'
        s += ');\n\n'
        return s

    def _edges(self):
        return 'edges\n(\n);\n\n'

    def _boundary(self):
        s = 'boundary\n(\n'
        # write out interfaces
        self.patch_pairs = []
        for i in range(self.Nlayers-1):
            name1 = 'interface{:d}{:d}'.format(i+1,i+2)
            name2 = 'interface{:d}{:d}'.format(i+2,i+1)
            self.patch_pairs.append((name1,name2))
            s += patch_def.format(name=name1,
                                  faceslist='            '+self.blocks[i].top())
            s += patch_def.format(name=name2,
                                  faceslist='            '+self.blocks[i].bottom())
        # lower boundary
        lowerfaces = '            '+self.blocks[0].lower()
        s += patch_def.format(name='lower',faceslist=lowerfaces)
        # upper boundary
        upperfaces = '            '+self.blocks[-1].upper()
        s += patch_def.format(name='upper',faceslist=upperfaces)
        # west boundary
        westfaces = '\n'.join([ '            '+self.blocks[i].west()
                                for i in range(self.Nlayers) ])
        s += patch_def.format(name='west',faceslist=westfaces)
        # east boundary
        eastfaces = '\n'.join([ '            '+self.blocks[i].east()
                                for i in range(self.Nlayers) ])
        s += patch_def.format(name='east',faceslist=eastfaces)
        # north boundary
        northfaces = '\n'.join([ '            '+self.blocks[i].north()
                                 for i in range(self.Nlayers) ])
        s += patch_def.format(name='north',faceslist=northfaces)
        # south boundary
        southfaces = '\n'.join([ '            '+self.blocks[i].south()
                                 for i in range(self.Nlayers) ])
        s += patch_def.format(name='south',faceslist=southfaces)
        s += ');\n\n'
        return s

    def _mergedPatchPairs(self):
        s = 'mergedPatchPairs\n(\n'
        for pair in self.patch_pairs:
            s += '    ' + str(pair) + '\n'
        s += ');\n'
        return s

    def write(self,fpath='blockMeshDict'):
        """Write out blockMeshDict, which should be placed into
        constant/polyMesh.
        """
        assert(len(self.Nx) == len(self.Ny) == len(self.Nz) ==
               len(self.z0) == len(self.z1) == len(self.simpleGradingZ))
        self.Nlayers = len(self.Nx)
        with open(fpath,'w') as f:
            f.write(blockMeshDict_header)
            f.write(self._vertices())
            f.write(self._blocks())
            f.write(self._edges())
            f.write(self._boundary())
            f.write(self._mergedPatchPairs())
            f.write(blockMeshDict_footer)
        print('Wrote '+fpath)


class hex(object):
    """Helper class for managing hex blocks"""

    def __init__(self, vertex0=0, N=(0,0,0), grading='simpleGrading', simpleGrading=(1,1,1)):
        self.vertices = np.arange(vertex0,vertex0+8)
        self.Ncells = N
        self.grading = grading
        self.simpleGrading = simpleGrading

    def write(self):
        return 'hex ({:s}) ({:s}) {:s} ({:s})'.format(
                ' '.join(['{:d}'.format(val) for val in self.vertices]),
                ' '.join(['{:d}'.format(val) for val in self.Ncells]),
                self.grading,
                ' '.join(['{:g}'.format(val) for val in self.simpleGrading])
                )

    def _facelist(self,v):
        return '('+' '.join([str(i) for i in v])+')'

    def west(self):
        return self._facelist(self.vertices[[0,4,7,3]])

    def east(self):
        return self._facelist(self.vertices[[1,2,6,5]])

    def south(self):
        return self._facelist(self.vertices[[0,1,5,4]])

    def north(self):
        return self._facelist(self.vertices[[2,3,7,6]])

    def lower(self):
        return self._facelist(self.vertices[[0,3,2,1]])

    def upper(self):
        return self._facelist(self.vertices[[4,5,6,7]])

