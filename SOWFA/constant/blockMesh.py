#
# Module for assisting with blockMesh grid generation
#
# Written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import numpy as np


"""General mesh functions"""

def grow_mesh(N, d0, r, verbose=True):
    """Calculate the cell spacings given a number of cells (N), an 
    initial spacing (d0), and a normal growth rate (r).

    Spacings are calculated as follows:
        d_i = d_0 * r^i,  for i in 0 .. N-1
    """
    d = d0 * np.array([r**i for i in range(N)])
    if verbose: print('spacings : {:g} .. {:g}'.format(d0,d[-1]))
    return d

def points_from(spacings, x0=0.0):
    x = x0 * np.ones((len(spacings)+1,))
    x[1:] += np.cumsum(spacings)
    return x

def spacings(L, N, ratio=1, verbose=True):
    """Calculate the blockmesh spacings that would result given a length
    (L), number of cells (N), and simpleGrading ratio.
    """
    r = np.exp(np.log(ratio) / (N-1))
    if verbose: print('growth rate = ',r)
    d0 = L / np.sum([r**i for i in range(N)])
    return grow_mesh(N, d0, r, verbose)

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
    # we assume that the estimated number of cells is pretty close to
    # what it should be and leave this fixed; since the growth rate
    # r = f(d0,d1,N), we should update it assuming the other parameters
    # are also fixed
    r_approx = np.exp(np.log(d/d0) / (N-1))
    print('actual growth rate: ',r_approx)
    # we're basically done at this point...
    L1 = np.sum(grow_mesh(N,d0,r_approx))
    print('expected layer height : ',L1)
    # but the resulting distance spanned by the points probably isn't a nice
    # floating point number or we guessed an incorrect layer height
    if L is None:
        L = np.round(L1/round_to) * round_to
    d0_approx = L / np.sum([r_approx**i for i in range(N)])
    print('adjusted initial spacing to obtain a height of {:g}: {:g}'.format(L,d0_approx))
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
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1.0;
"""

patch_def = """    {name:s}
    {{
        type {ptype:s};
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

        Note: zMax is only used for generating a single-block uniform
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

        # blocks (i.e., layers, only hexes for now)
        self.vertex0 = 0
        self.blocks = []
        self.patchpairs = []

        # descriptors of each layer
        self.Nx = []
        self.Ny = []
        self.Nz = []
        self.z0 = []
        self.z1 = []
        self.simpleGradingZ = []

    def size(self):
        return np.sum([ block.size() for block in self.blocks ])

    def __repr__(self):
        Nlayers = len(self.Nx)
        zMin = self.zMin
        zMax = self.zMax
        if Nlayers > 0:
            zMin = self.z0[0]
            zMax = self.z1[-1]
        s = 'Mesh bounding box corners: ({:g} {:g} {:g}) ({:g} {:g} {:g})\n'.format(
                self.xMin, self.yMin, zMin,
                self.xMax, self.yMax, zMax)
        if Nlayers == 0:
            s += '  no layers added; call generate_uniform_grid or generate_layer with add=True'
        for i in range(Nlayers):
            dx = (self.xMax - self.xMin) / self.Nx[i]
            dy = (self.yMax - self.yMin) / self.Ny[i]
            Lz = self.z1[i] - self.z0[i]
            if self.simpleGradingZ[i] == 1:
                dz = Lz / self.Nz[i]
                dzstr = '{:g}'.format(dz)
            else:
                d = spacings(Lz, self.Nz[i], self.simpleGradingZ[i], verbose=False)
                dzstr = '{:g}..{:g}'.format(d[0],d[-1])
            s+= '  layer {:d} : z=({:g} {:g}) N=({:d} {:d} {:d}) spacings=({:g} {:g} {:s})\n'.format(
                    i, self.z0[i], self.z1[i],
                    self.Nx[i], self.Ny[i], self.Nz[i],
                    dx, dy, dzstr)
        if Nlayers > 0:
            s += 'Total cells: {:d}\n'.format(self.size())
        return s
            

    def generate_uniform_grid(self):
        """Wrapper around generate_layer() for a simple single-block
        mesh with uniform grid spacing throughout
        """
        self.generate_layer(self.zMin, self.zMax, self.dx0, self.dy0, self.dz0,
                             add=True)

    def generate_layer(self,z0,z1,dx,dy,dz0,dz1=None,r=1,add=False):
        """Generate new layer between height z0 and z1, with grid
        spacings dx, dy, and dz.

        Optionally, a growth rate r > 1 may be specified, in which case
        the target final grid spacing (dz1) should be specified.

        The new layer is added to the domain until the add=True is
        specified.
        """
        Nx = int(np.ceil((self.xMax-self.xMin)/dx))
        x0 = self.xMin
        x1 = x0 + Nx*dx
        if x1 > self.xMax:
            self.xMax = x1
            print('updated xMax to',self.xMax)

        Ny = int(np.ceil((self.yMax-self.yMin)/dy))
        y0 = self.yMin
        y1 = y0 + Ny*dy
        if y1 > self.yMax:
            self.yMax = y1
            print('updated yMax to',self.yMax)

        if r==1:
            Nz = int(np.ceil((z1-z0)/dz0))
            new_z1 = z0 + Nz*dz0
            if new_z1 > z1:
                z1 = new_z1
                print('updated top of layer to',z1)
            ratio = 1.0
            print('constant vertical spacing layer, Nz=',Nz)
        else:
            assert(dz1 is not None)
            Lz,Nz,ratio = estimate(dz1,dz0,r,L=(z1-z0))
            assert(Lz == z1-z0)
        
        if add:
            # done designing; want to add this layer
            if (len(self.Nx) == 0) or (z0 == self.z1[-1]): 
                # this is either the first layer or it stacks nicely on top of
                # the previous layer
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

                print('Added layer',len(self.Nx))
            else:
                print('Invalid heights, layer not added')

    def _constants(self):
        s = '\n'
        for const in ['xMin','xMax','yMin','yMax']:
            s += '{:s}\t\t{:g};\n'.format(const,getattr(self,const))
        for i in range(self.Nlayers):
            s += 'zMin{:d}\t\t{:g};\n'.format(i,self.z0[i])
            s += 'zMax{:d}\t\t{:g};\n'.format(i,self.z1[i])
        s += '\n'
        return s

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
        self.patchpairs = []
        for i in range(self.Nlayers-1):
            name1 = 'interface{:d}{:d}'.format(i+1,i+2)
            name2 = 'interface{:d}{:d}'.format(i+2,i+1)
            self.patchpairs.append((name1,name2))
            s += patch_def.format(name=name1, ptype='patch',
                                  faceslist='            '+self.blocks[i].upper())
            s += patch_def.format(name=name2, ptype='patch',
                                  faceslist='            '+self.blocks[i+1].lower())
        # lower boundary
        lowerfaces = '            '+self.blocks[0].lower()
        s += patch_def.format(name='lower',ptype='wall',faceslist=lowerfaces)
        # upper boundary
        upperfaces = '            '+self.blocks[-1].upper()
        s += patch_def.format(name='upper',ptype='patch',faceslist=upperfaces)
        # west boundary
        westfaces = '\n'.join([ '            '+self.blocks[i].west()
                                for i in range(self.Nlayers) ])
        s += patch_def.format(name='west',ptype='patch',faceslist=westfaces)
        # east boundary
        eastfaces = '\n'.join([ '            '+self.blocks[i].east()
                                for i in range(self.Nlayers) ])
        s += patch_def.format(name='east',ptype='patch',faceslist=eastfaces)
        # north boundary
        northfaces = '\n'.join([ '            '+self.blocks[i].north()
                                 for i in range(self.Nlayers) ])
        s += patch_def.format(name='north',ptype='patch',faceslist=northfaces)
        # south boundary
        southfaces = '\n'.join([ '            '+self.blocks[i].south()
                                 for i in range(self.Nlayers) ])
        s += patch_def.format(name='south',ptype='patch',faceslist=southfaces)
        s += ');\n\n'
        return s

    def _mergePatchPairs(self):
        s = 'mergePatchPairs\n(\n'
        for pair in self.patchpairs:
            s += '    ({:s} {:s})\n'.format(*pair)
        s += ');\n'
        return s

    def write(self,fpath='blockMeshDict'):
        """Write out blockMeshDict, which should be placed into
        constant/polyMesh.
        """
        if len(self.blocks) == 0:
            raise IndexError('No blocks have been generated')
        assert(len(self.Nx) == len(self.Ny) == len(self.Nz) ==
               len(self.z0) == len(self.z1) == len(self.simpleGradingZ))
        self.Nlayers = len(self.Nx)
        with open(fpath,'w') as f:
            f.write(blockMeshDict_header)
            f.write('\n\n/*\n' +
                    ' * Generated using NWTC/datatools python library\n' +
                    ' * https://github.com/NWTC/datatools/blob/master/SOWFA/constant/blockMesh.py\n' +
                    ' * \n\n' +
                    self.__repr__() +
                    '*/\n\n')
            f.write(self._constants())
            f.write(self._vertices())
            f.write(self._blocks())
            f.write(self._edges())
            f.write(self._boundary())
            f.write(self._mergePatchPairs())
            f.write(blockMeshDict_footer)
        print('Wrote '+fpath)


class hex(object):
    """Helper class for managing hex blocks"""

    def __init__(self, vertex0=0, N=(0,0,0), grading='simpleGrading', simpleGrading=(1,1,1)):
        self.vertices = np.arange(vertex0,vertex0+8)
        self.Ncells = N
        self.grading = grading
        self.simpleGrading = simpleGrading

    def size(self):
        return np.prod(self.Ncells)

    def write(self):
        return 'hex ({:s}) ({:s}) {:s} ({:s})'.format(
                ' '.join(['{:d}'.format(val) for val in self.vertices]),
                ' '.join(['{:d}'.format(val) for val in self.Ncells]),
                self.grading,
                ' '.join(['{:.12g}'.format(val) for val in self.simpleGrading])
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

