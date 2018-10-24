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
    source_types = ['background_box','turbine_box','turbine_cylinder']

    def __init__(self,sources=[],perturb=0.01,rotation=0.0,
                 upstream=5.0,downstream=10.0,width=3.0,height=3.0,
                 streamwise_buffer=1.0,lateral_buffer=1.0,vertical_buffer=1.0):
        """Object for generating a series of topoSetDict files (e.g.,
        for refinement). Level 1 is the finest level, and successive
        refinements should be performed sequentially starting from the
        coarsest level down to 1.

        Note that defaults are only specified for turbine refinement
        regions; the actual dimensional values must be specified for 
        each background refinement region.
        
        Inputs
        ------
        sources : list of str
            Describes the cellSet sources (either "box" or "cylinder")
            used in describing the refinement regions, the length of
            which corresponds to the number of refinement levels.
        perturb : float, optional
            A perturbation of the refinement boxes to keep the
            boundaries off of cell faces. [m]

        Default turbine box refinement parameters
        -----------------------------------------
        rotation : float, optional
            Angle about the z-axis to rotate the refinement region (NOT
            the compass wind direction); if None, then mean_rotation is
            used. [deg]
        upstream : float, optional
            Distance (in diameters) upstream of the turbine where the
            inner refinement box starts.
        downstream : float, optional
            Distance (in diameters) downstream of the turbine where the
            inner refinement box ends.
        width : float, optional
            The overall width (in diameters) of the inner refinement
            box; need to account for horizontal wake meandering
            downstream.
        height : float, optional
            The overall height (in diameters) of the inner refinement
            box; need to account for vertical wake meandering
            downstream.
        streamwise_buffer : float, optional
            Size of buffer region (in diameters) between refinement
            levels in the upstream/downstream directions
        lateral_buffer : float, optional
            Size of buffer region (in diameters) between refinement
            levels in the lateral directions
        vertical_buffer : float, optional
            Size of buffer region (in diameters) between refinement
            levels in the vertical directions
        """
        self.sources = sources
        self.Nlevels = len(sources)
        self._check_specified_sources()

        # set defaults
        self.refinement = dict(
            rotation=rotation*np.pi/180.0,
            upstream=upstream, downstream=downstream,
            width=width, height=height,
            streamwise_buffer=streamwise_buffer,
            lateral_buffer=lateral_buffer,
            vertical_buffer=vertical_buffer
        )

        # definitions for each background region
        self.LLcorner = []
        self.bkg_rotation = []
        self.bkg_length = []
        self.bkg_width = []
        self.bkg_height = []
        self.bkg_xbuffer = []
        self.bkg_ybuffer = []
        self.bkg_zbuffer = []

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
            

    def __repr__(self):
        s = '{:d} refinement levels : {:s}'.format(self.Nlevels,
                                                   str(self.sources))
        for iturb,loc in enumerate(self.base_location):
            s += '\nturbine {:d} at {:s} rotated {:g} deg'.format(
                    iturb+1, str(loc), 180./np.pi*self.rotation[iturb])
        return s


    def _check_specified_sources(self):
        defined = [ (sourcename in self.source_types)
                    for sourcename in self.sources ]
        assert(all(defined))

    def add_background(self, LLcorner=(0,0,0),
            length=0.0, width=0.0, height=0.0, rotation=0.0,
            streamwise_buffer=50.0, lateral_buffer=50.0, vertical_buffer=50.0):
        """Add refinement box at location specified by the lower-left
        corner 'LLcorner' with dimensions given by Lx, Ly, Lz.
        
        By default, the box is aligned with the x-, y-, and z- axes.

        Background box parameters
        -------------------------
        LLcorner : array-like
            xyz coordinates of lower-left corner of the refinement box
            [m]
        length : float
            Length of box in the x-direction. [m]
        width : float
            Width of box in the y-direction. [m]
        height : float
            Height of box in the x-direction. [m]
        rotation : float, optional
            Angle about the z-axis to rotate the refinement region. [deg]
        streamwise_buffer : float, optional
            Size of buffer region between refinement levels in the upstream
            and downstream directions. [m]
        lateral_buffer : float, optional
            Size of buffer region between refinement levels in the lateral
            directions. [m]
        vertical_buffer : float, optional
            Size of buffer region between refinement levels in the vertical
            direction. [m]
        """
        assert((length > 0) and (width > 0) and (height > 0))
        self.LLcorner.append(LLcorner)
        if 'rotation' is None:
            rotation = self.refinement['rotation']
        self.bkg_rotation.append(rotation)
        self.bkg_length.append(length)
        self.bkg_width.append(width)
        self.bkg_height.append(height)
        self.bkg_xbuffer.append(streamwise_buffer)
        self.bkg_ybuffer.append(lateral_buffer)
        self.bkg_zbuffer.append(vertical_buffer)

    def add_turbine(self, base_location=(1000,1000,0), D=126.0, **kwargs):
        """Add turbine at specified 'base_location' with diameter 'D'.

        Note that each turbine is associated with a set of topoSet
        refinement sources.

        Turbine parameters
        ------------------
        base_location : array-like
            xyz coordinates of turbine base. [m]
        D : float
            Rotor dameter, used as a reference length. [m]
        **kwargs :
            Optional refinement parameters to override defaults
            specified during initialization.

        """
        self.base_location.append(base_location)
        self.diameter.append(D)
        if 'rotation' in kwargs:
            ang = np.pi/180. * kwargs['rotation']
        else:
            ang = self.refinement['rotation']
        def get_param(param):
            return kwargs.get(param, self.refinement[param])
        self.rotation.append(ang)
        self.upstream.append(get_param('upstream'))
        self.downstream.append(get_param('downstream'))
        self.width.append(get_param('width'))
        self.height.append(get_param('height'))
        self.xbuffer.append(get_param('streamwise_buffer'))
        self.ybuffer.append(get_param('lateral_buffer'))
        self.zbuffer.append(get_param('vertical_buffer'))


    def write(self,prefix='topoSetDict.local'):
        for ilevel in range(self.Nlevels):
            fname = '{:s}.{:d}'.format(prefix,ilevel+1)
            sourcename = self.sources[ilevel]
            source = getattr(self,'_write_'+sourcename)
#            print('Writing {:s} dict : {:s}'.format(sourcename,fname))
            # Get the effective level; if sources==['cylinder','box','box'],
            # then the zero-indexed ilevel==1 (corresponds to an overall
            # refinement level of 2, and a box-refinement level of 1)
            efflevel = ilevel
            for i in range(ilevel):
                if not self.sources[i] == self.sources[ilevel]:
                    efflevel -= 1
            print('Writing {:s} dict, level {:d} : {:s}'.format(sourcename,
                                                                efflevel,
                                                                fname))
            # write out topoSetDict.*
            with open(fname,'w') as f:
                f.write(topoSetDict_header.format(fname=fname))
                if sourcename.startswith('turbine'):
                    for iturb in range(len(self.base_location)):
                        f.write(source(iturb,efflevel)) # levels are 0-indexed
                elif sourcename.startswith('background'):
                    for ibkg in range(len(self.LLcorner)):
                        f.write(source(ibkg,efflevel)) # levels are 0-indexed
                else:
                    print("Unrecognized source type for source '{:s}'".format(sourcename))
                f.write(topoSetDict_footer)


    def _write_background_box(self,ibkg,ilevel):
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
        if ibkg == 0:
            action = 'new'
        else:
            action = 'add'
        length = self.bkg_length[ibkg]
        width = self.bkg_width[ibkg]
        height = self.bkg_height[ibkg]
        xbuff = self.bkg_xbuffer[ibkg]
        ybuff = self.bkg_ybuffer[ibkg]
        zbuff = self.bkg_zbuffer[ibkg]
        LLcorner = self.LLcorner[ibkg]
        ang = self.bkg_rotation[ibkg]
        # origin (x,y,z)
        x0 = -ilevel*xbuff
        y0 = -ilevel*ybuff
        x = x0 * np.cos(ang) - y0 * np.sin(ang) + LLcorner[0]
        y = x0 * np.sin(ang) + y0 * np.cos(ang) + LLcorner[1]
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


    def _write_turbine_box(self,iturb,ilevel):
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


    def _write_turbine_cylinder(self,iturb,ilevel):
        template = """    {{
        name         local;
        type         cellSet;
        action       {action:s};
        source       ;
        sourceInfo
        {{
        }}
    }}
"""
        if iturb == 0:
            action = 'new'
        else:
            action = 'add'
        raise NotImplementedError('_write_cylinder')

