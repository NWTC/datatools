#
# Module to generate topoSet input dictionaries
#
# Written by Eliot Quon (eliot.quon@nrel.gov)
#
import numpy as np

topoSetDict_header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.0                                   |
|   \\\\  /    A nd           | Web:      http://www.OpenFOAM.org               |
|    \\\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      {fname};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

actions
(
""" 
topoSetDict_footer = """);

// ************************************************************************* //"""


class TopoSetDict(object):
    source_types = ['box','cylinder']

    def __init__(self,refinement=[],mean_rotation=0.0,perturb=0.01):
        """Object for generating a series of topoSetDict files for
        refinement. Level 1 is the finest level, and successive
        refinements should be performed sequentially starting from the
        coarsest level down to 1.
        
        Inputs
        ------
        refinement : list of str
            Describes the cellSet sources (either "box" or "cylinder")
            used in describing the refinement regions, the length of
            which corresponds to the number of refinement levels.
        mean_rotation : float, optional
            Default rotation of refinement regions about z-axis using
            the right-hand rule (NOT compass wind direction). [deg]
        perturb : float, optional
            A perturbation of the refinement boxes to keep the
            boundaries off of cell faces. [m]
        """
        self.refinement = refinement
        self.Nlevels = len(refinement)
        self._check_specified_refinement()
        self.rotation0 = mean_rotation * np.pi/180.0

        # definitions for each turbine
        self.base_location = []
        self.rotation = []
        self.diameter = []
        self.upstream = []
        self.downstream = []
        self.width = []
        self.height = []
        self.xbuffer = []
        self.ybuffer = []
        self.zbuffer = []
            

    def _check_specified_refinement(self):
        defined = [ (refinetype in self.source_types)
                    for refinetype in self.refinement ]
        assert(all(defined))


    def __repr__(self):
        s = '{:d} refinement levels : {:s}'.format(self.Nlevels,
                                                   str(self.refinement))
        for iturb,loc in enumerate(self.base_location):
            s += '\nturbine {:d} at {:s} rotated {:g} deg'.format(
                    iturb+1, str(loc), 180./np.pi*self.rotation[iturb])
        return s


    def add_turbine(self, base_location=(1000,1000,0),
                rotation=None,
                D=126.0,
                upstream=5.0,
                downstream=10.0,
                width=3.0,
                height=3.0,
                streamwise_buffer=1.0,
                lateral_buffer=1.0,
                vertical_buffer=1.0):
        """Add turbine at specified 'base_location' with diameter 'D'.
        The 'upstream', 'downstream', 'width', and 'height' lists should
        all have the same length as the self.refinement list.

        Refinement parameters
        ---------------------
        base_location : array-like
            xyz coordinates of turbine base. [m]
        rotation : float, optional
            Angle about the z-axis to rotate the refinement region (NOT
            the compass wind direction); if None, then mean_rotation is
            used. [deg]
        D : float
            Rotor dameter, used as a reference length. [m]
        upstream : float
            Distance (in diameters) upstream of the turbine where the
            inner refinement box starts.
        downstream : float
            Distance (in diameters) downstream of the turbine where the
            inner refinement box ends.
        width : float
            The overall width (in diameters) of the inner refinement
            box; need to account for horizontal wake meandering
            downstream.
        height : float
            The overall height (in diameters) of the inner refinement
            box; need to account for vertical wake meandering
            downstream.
        streamwise_buffer : float
            Size of buffer region (in diameters) between refinement
            levels in the upstream/downstream directions
        lateral_buffer : float
            Size of buffer region (in diameters) between refinement
            levels in the lateral directions
        vertical_buffer : float
            Size of buffer region (in diameters) between refinement
            levels in the vertical directions
        """
        self.base_location.append(base_location)
        if rotation is None:
            ang = self.rotation0
        else:
            ang = np.pi/180. * rotation
        self.rotation.append(ang)            
        self.diameter.append(D)
        self.upstream.append(upstream)
        self.downstream.append(downstream)
        self.width.append(width)
        self.height.append(height)
        self.xbuffer.append(streamwise_buffer)
        self.ybuffer.append(lateral_buffer)
        self.zbuffer.append(vertical_buffer)


    def write(self,prefix='topoSetDict.local'):
        for ilevel in range(self.Nlevels):
            fname = '{:s}.{:d}'.format(prefix,ilevel+1)
            refinetype = self.refinement[ilevel]
            source = getattr(self,'_write_'+refinetype)
            print('Writing {:s} dict : {:s}'.format(refinetype,fname))
            # Get the effective level; if refinement==['cylinder','box','box'],
            # then the zero-indexed ilevel==1 (corresponds to an overall
            # refinement level of 2, and a box-refinement level of 1)
#            efflevel = ilevel
#            for i in range(ilevel):
#                if not self.refinement[i] == self.refinement[ilevel]:
#                    efflevel -= 1
#            print('Writing {:s} dict, level {:d} : {:s}'.format(refinetype,
#                                                                efflevel,
#                                                                fname))
            # write out topoSetDict.*
            with open(fname,'w') as f:
                f.write(topoSetDict_header.format(fname=fname))
                for iturb in range(len(self.base_location)):
                    #f.write(source(iturb,efflevel))
                    f.write(source(iturb,ilevel)) # ilevel is 0-indexed
                f.write(topoSetDict_footer)


    def _write_box(self,iturb,ilevel):
        template = """    {{
        name         local;
        type         cellSet;
        action       {action:s};
        source       rotatedBoxToCell;
        sourceInfo
        {{
            origin ( {x0:g} {y0:g} {z0:g} );
            i      ( {ix:g} {iy:g} {iz:g} );
            j      ( {jx:g} {jy:g} {jz:g} );
            k      ( {kx:g} {ky:g} {kz:g} );
        }}
    }}
"""
        if iturb == 0:
            action = 'new'
        else:
            action = 'add'
        Lref = self.diameter[iturb]
        upstream = self.upstream[iturb] * Lref
        downstream = self.downstream[iturb] * Lref
        length = upstream + downstream
        width = self.width[iturb] * Lref
        height = self.height[iturb] * Lref
        xbuff = self.xbuffer[iturb] * Lref
        ybuff = self.ybuffer[iturb] * Lref
        zbuff = self.zbuffer[iturb] * Lref
        base = self.base_location[iturb]
        ang = self.rotation[iturb]
        # origin (x,y,z)
        x0 = -upstream - ilevel*xbuff
        y0 = -0.5*width - ilevel*ybuff
        x = x0 * np.cos(ang) - y0 * np.sin(ang) + base[0]
        y = x0 * np.sin(ang) + y0 * np.cos(ang) + base[1]
        z = base[2]
        # box dimensions
        ix = (length + 2*ilevel*xbuff) * np.cos(ang)
        iy = (length + 2*ilevel*xbuff) * np.sin(ang)
        iz = 0.0
        jx = -(width + 2*ilevel*ybuff) * np.sin(ang)
        jy =  (width + 2*ilevel*ybuff) * np.cos(ang)
        jz = 0.0
        kx = 0.0
        ky = 0.0
        kz = height + ilevel*zbuff
        return template.format(action=action,
                x0=x, y0=y, z0=z,
                ix=ix, iy=iy, iz=iz,
                jx=jx, jy=jy, jz=jz,
                kx=kx, ky=ky, kz=kz)


    def _write_cylinder(self,f,iturb,loc):
        template = """    {{
        name         local;
        type         cellSet;
        action       {action:s};
        source       rotatedBoxToCell;
        sourceInfo
        {{
            origin ( {x0:g} {y0:g} {z0:g} );
            i      ( {ix:g} {iy:g} {iz:g} );
            j      ( {jx:g} {jy:g} {jz:g} );
            k      ( {kx:g} {ky:g} {kz:g} );
        }}
    }}
"""
        if iturb == 0:
            action = 'new'
        else:
            action = 'add'
        raise NotImplementedError('_write_cylinder')

