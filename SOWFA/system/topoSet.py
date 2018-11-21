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
    defaults = dict(
        rotation=0.0,
        upstream=5.0, downstream=10.0,
        width=3.0, height=3.0,
        xbuffer=1.0, ybuffer=1.0, zbuffer=1.0,
        rbuffer=0.5
    )

    def __init__(self,sources=[],perturb=0.01,**kwargs):
        """Object for generating a series of topoSetDict files (e.g.,
        for refinement). Level 1 is the finest level, and successive
        refinements should be performed sequentially starting from the
        coarsest level down to 1.

        Note that defaults are only specified for turbine refinement
        regions; the actual dimensional values must be specified for 
        each background refinement region.

        Optional keyword arguments are passed to setup()
        
        Inputs
        ------
        sources : list of str
            cellSet sources (available in sources in self.source_types)
            that describe the refinement regions, the length of which
            corresponds to the number of refinement levels.
        perturb : float, optional
            A perturbation of the refinement boxes to keep the
            boundaries off of cell faces. [m]
        """
        self.sources = sources
        self._check_specified_sources()

        self.refinement = dict()
        if len(kwargs) == 0:
            print('using defaults')
        self.setup(**kwargs)

        # definitions for each background region
        self.bkg_LLcorner = []
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
        self.zhub = []
        self.upstream = []
        self.downstream = []
        self.width = []
        self.height = []
        self.xbuffer = []
        self.ybuffer = []
        self.zbuffer = []
        self.rbuffer = []
            

    def __repr__(self):
        s = '{:d} refinement levels : {:s}'.format(len(self.sources),
                                                   str(self.sources))
        for ibkg,loc in enumerate(self.bkg_LLcorner):
            s += '\nbackground region {:d} at {:s} rotated {:g} deg'.format(
                    ibkg+1, str(loc), self.bkg_rotation[ibkg])
        for iturb,loc in enumerate(self.base_location):
            s += '\nturbine {:d} at {:s} rotated {:g} deg'.format(
                    iturb+1, str(loc), self.rotation[iturb])
        return s


    def _check_specified_sources(self):
        """Check that all specified sources are defined, i.e., they have
        associated _write_* functions defined.
        """
        defined = [ (sourcename in self.source_types)
                    for sourcename in self.sources ]
        assert(all(defined))

    
    def setup(self,**kwargs):
        """
        Refinement Parameters
        ---------------------
        rotation : float
            Angle about the z-axis to rotate the refinement region (NOT
            the compass wind direction); if None, then mean_rotation is
            used. [deg]
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
        xbuffer : float
            Size of buffer region (in diameters) between refinement
            levels in the upstream/downstream directions
        ybuffer : float
            Size of buffer region (in diameters) between refinement
            levels in the lateral directions
        zbuffer : float
            Size of buffer region (in diameters) between refinement
            levels in the vertical directions
        rbuffer : float
            Used to set size of cylindrical refinement region; cylinder
            level i has diameter (1 + (i+1)*rbuffer)*D for
            i=0,1,...
        """
        for key,defval in self.defaults.items():
            val = kwargs.get(key, defval)
            #if not hasattr(val,'__iter__')
            self.refinement[key] = val
        print('refinement parameters:',self.refinement)


    def add_background_box(self, LLcorner=(0,0,0),
            length=0.0, width=0.0, height=0.0, rotation=0.0,
            xbuffer=50.0, ybuffer=50.0, zbuffer=50.0):
        """Add refinement box at location specified by the lower-left
        corner 'LLcorner' with dimensions given by Lx, Ly, Lz.

        Background refinement boxes are only implemented for a single level.
        
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
        xbuffer : float, optional
            Size of buffer region between refinement levels in the upstream
            and downstream directions. [m]
        ybuffer : float, optional
            Size of buffer region between refinement levels in the lateral
            directions. [m]
        zbuffer : float, optional
            Size of buffer region between refinement levels in the vertical
            direction. [m]
        """
        if not self.sources[-1] == 'background_box':
            self.sources.append('background_box')
        assert((length > 0) and (width > 0) and (height > 0))
        self.bkg_LLcorner.append(LLcorner)
        if 'rotation' is None:
            rotation = self.refinement['rotation']
        self.bkg_rotation.append(rotation)
        self.bkg_length.append(length)
        self.bkg_width.append(width)
        self.bkg_height.append(height)
        self.bkg_xbuffer.append(xbuffer)
        self.bkg_ybuffer.append(ybuffer)
        self.bkg_zbuffer.append(zbuffer)

    def add_turbine(self,
            base_location=(1000,1000,0), D=126.0, zhub=None,
            **kwargs):
        """Add turbine at specified 'base_location' with diameter 'D'
        and orientation 'rotation'.

        Note that each turbine is associated with a set of topoSet
        refinement sources.

        Turbine parameters
        ------------------
        base_location : array-like
            xyz coordinates of turbine base. [m]
        D : float
            Rotor diameter, used as a reference length. [m]
        zhub : float, optional
            Hub height, used for cylindrical refinement; if None, set
            to rotor diameter. [m]
        **kwargs :
            Optional refinement parameters to override defaults
            specified during initialization.
        """
        self.base_location.append(base_location)
        self.diameter.append(D)
        if zhub is None:
            zhub = D
        self.zhub.append(zhub)
        if 'rotation' in kwargs:
            ang_deg = kwargs['rotation']
        else:
            ang_deg = self.refinement['rotation']
        def get_param(param):
            return kwargs.get(param, self.refinement[param])
        self.rotation.append(np.radians(ang_deg))
        self.upstream.append(get_param('upstream'))
        self.downstream.append(get_param('downstream'))
        self.width.append(get_param('width'))
        self.height.append(get_param('height'))
        self.xbuffer.append(get_param('xbuffer'))
        self.ybuffer.append(get_param('ybuffer'))
        self.zbuffer.append(get_param('zbuffer'))
        self.rbuffer.append(get_param('rbuffer'))


    def estimate_mesh_size(self,initial_size=None,ds0=10.0):
        """NOTE: this is experimental and not very accurate!"""
        if initial_size is None:
            raise ValueError('specify initial cell count or list of dimensions')
        if hasattr(initial_size,'__iter__'):
            initial_size = int(np.prod(initial_size))
            print('calculated initial cell count: {:d}'.format(initial_size))
        Nlevels = len(self.sources)
        vol = np.zeros(Nlevels)
        for ilevel in range(Nlevels):
            sourcename = self.sources[ilevel]
            efflevel = ilevel
            for i in range(ilevel):
                if not self.sources[i] == self.sources[ilevel]:
                    efflevel -= 1
            if sourcename == 'background_box':
                print('{:d}: background regions {:d}'.format(ilevel,efflevel))
                for ibkg in range(len(self.bkg_LLcorner)):
                    length = self.bkg_length[ibkg]
                    width = self.bkg_width[ibkg]
                    height = self.bkg_height[ibkg]
                    xbuff = self.bkg_xbuffer[ibkg]
                    ybuff = self.bkg_ybuffer[ibkg]
                    zbuff = self.bkg_zbuffer[ibkg]
                    length += 2*efflevel*xbuff
                    width += 2*efflevel*ybuff
                    height += efflevel*zbuff
                    vol[ilevel] += length*width*height
                    #print('  box{:d} vol  = {:g}'.format(ibkg,length*width*height))
            elif sourcename == 'box':
                print('{:d}: turbine box regions {:d}'.format(ilevel,efflevel))
                for iturb in range(len(self.base_location)):
                    Lref = self.diameter[iturb]
                    upstream = self.upstream[iturb] * Lref
                    downstream = self.downstream[iturb] * Lref
                    length = upstream + downstream
                    width = self.width[iturb] * Lref
                    height = self.height[iturb] * Lref
                    xbuff = self.xbuffer[iturb] * Lref
                    ybuff = self.ybuffer[iturb] * Lref
                    zbuff = self.zbuffer[iturb] * Lref
                    length += 2*efflevel*xbuff
                    width += 2*efflevel*ybuff
                    height += efflevel*zbuff
                    vol[ilevel] += length*width*height
            elif sourcename == 'cylinder':
                print('{:d}: turbine cylinder regions {:d}'.format(ilevel,efflevel))
                for iturb in range(len(self.base_location)):
                    Lref = self.diameter[iturb]
                    R = (1.0 + (efflevel+1)*self.rbuffer[iturb]) * Lref/2
                    upstream = self.upstream[iturb] * Lref
                    downstream = self.downstream[iturb] * Lref
                    xbuff = self.xbuffer[iturb] * Lref
                    length = upstream + downstream + 2*efflevel*xbuff
                    vol[ilevel] += length * np.pi*R**2
        lastcount = initial_size
        ds = float(ds0)
        for ilevel in range(Nlevels):
            print('volume, ds : {:g} {:f}'.format(vol[ilevel],ds))
            approx_cells = int(vol[ilevel] / ds**3)
            print('approx number of selected cells : {:d}'.format(approx_cells))
            newcount = lastcount + 7*approx_cells
            print('refinement level {:d} : increase from {:d} to {:d}'.format(
                    ilevel+1, lastcount, newcount))
            lastcount = newcount
            ds /= 2


    def write(self,prefix='topoSetDict.local'):
        for ilevel in range(len(self.sources)):
            fname = '{:s}.{:d}'.format(prefix,ilevel+1)
            sourcename = self.sources[ilevel]
            source = getattr(self,'_write_'+sourcename)
            # Get the effective level; e.g., if
            #   sources==['cylinder','box','box'],
            # then the zero-indexed ilevel==1 (corresponds to an overall
            # refinement level of 2, and a box-refinement level of 1)
            efflevel = ilevel
            for i in range(ilevel):
                if not self.sources[i] == self.sources[ilevel]:
                    efflevel -= 1
            # write out topoSetDict.*
            with open(fname,'w') as f:
                f.write(topoSetDict_header.format(fname=fname))
                if sourcename.startswith('background'):
                    print('{:d}: Writing {:s} dict, level {:d} to {:s}'.format(
                            ilevel, sourcename, efflevel, fname))
                    for ibkg in range(len(self.bkg_LLcorner)):
                        f.write(source(ibkg,efflevel)) # levels are 0-indexed
                elif sourcename in self.source_types:
                    print('{:d}: Writing {:s} dict, level {:d} to {:s}'.format(
                            ilevel, sourcename, efflevel, fname))
                    for iturb in range(len(self.base_location)):
                        f.write(source(iturb,efflevel)) # levels are 0-indexed
                else:
                    print("Unrecognized source type for source '{:s}'".format(sourcename))
                f.write(topoSetDict_footer)


    def _write_background_box(self,ibkg,ilevel):
        # Depends on bkg_length, bkg_width, bkg_height, bkg_xbuffer,
        # bkg_ybuffer, and bkg_zbuffer.
        template = """    {{
        name         local;
        type         cellSet;
        action       {action:s};
        source       rotatedBoxToCell;
        sourceInfo
        {{
            origin ( {x0:f} {y0:f} {z0:f} );
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
        LLcorner = self.bkg_LLcorner[ibkg]
        ang = self.bkg_rotation[ibkg]
        # origin (x,y,z)
        x0 = -ilevel*xbuff
        y0 = -ilevel*ybuff
        x = x0 * np.cos(ang) - y0 * np.sin(ang) + LLcorner[0]
        y = x0 * np.sin(ang) + y0 * np.cos(ang) + LLcorner[1]
        z = LLcorner[2]
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


    def _write_box(self,iturb,ilevel):
        # Write out topoSetDict for a turbine refinement box.
        # Depends on D, upstream, downstream, width, height,
        # streamwise_buffer, lateral_buffer, and vertical_buffer.
        template = """    {{
        name         local;
        type         cellSet;
        action       {action:s};
        source       rotatedBoxToCell;
        sourceInfo
        {{
            origin ( {x0:f} {y0:f} {z0:f} );
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


    def _write_cylinder(self,iturb,ilevel):
        # Write out topoSetDict for a turbine refinement cylinder.
        # Depends on D, radial_buffer, upstream, downstream, and
        # streamwise buffer.
        template = """    {{
        name         local;
        type         cellSet;
        action       {action:s};
        source       cylinderToCell;
        sourceInfo
        {{
            p1      ( {x1:f} {y1:f} {z1:f} );
            p2      ( {x2:f} {y2:f} {z2:f} );
            radius  {R:g};
        }}
    }}
"""
        if iturb == 0:
            action = 'new'
        else:
            action = 'add'
        Lref = self.diameter[iturb]
        R = (1.0 + (ilevel+1)*self.rbuffer[iturb]) * Lref/2
        upstream = self.upstream[iturb] * Lref
        downstream = self.downstream[iturb] * Lref
        xbuff = self.xbuffer[iturb] * Lref
        base = self.base_location[iturb]
        ang = self.rotation[iturb]
        # upstream point
        xu = -upstream - ilevel*xbuff
        x1 = xu * np.cos(ang) + base[0]
        y1 = xu * np.sin(ang) + base[1]
        xd = downstream + ilevel*xbuff
        x2 = xd * np.cos(ang) + base[0]
        y2 = xd * np.sin(ang) + base[1]
        z = self.zhub[iturb] + base[2]
        return template.format(action=action,
                x1=x1, y1=y1, z1=z,
                x2=x2, y2=y2, z2=z,
                R=R)


    def _write_NEW(self,i,ilevel):
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

