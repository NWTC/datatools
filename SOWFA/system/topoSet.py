#
# Module to generate topoSet input dictionaries
#
# Written by Eliot Quon (eliot.quon@nrel.gov)
#
import os,sys
import numpy as np
import matplotlib.pyplot as plt

def plot(fpath,plane='xy',verbose=False,IDs=None,**plot_kwargs):
    """Plot an existing topoSetDict (only tested for topoSetDict files
    generated using the TopoSetDict class. One-index list of IDs may be
    provided to plot only specific turbine refinement regions.
    """
    with open(fpath,'r') as f:
        while not f.readline().strip() == 'actions': continue
        f.readline()
        def read_block():
            block = ''
            inputs = {}
            while (block=='') or \
                    (not block.count('{') == block.count('}')):
                line = f.readline()
                if line == '':
                    # EOF, but probably shouldn't end up here
                    return None
                line = line.strip()
                if line.startswith('//'):
                    continue
                block += line
                if line.endswith(';'):
                    # parse it
                    line = line.rstrip(';')
                    if '(' in line:
                        assert(line.count('(') == line.count(')'))
                        line = line.replace('(','').replace(')','')
                        is_array = True
                    else:
                        is_array = False
                    line = line.split()
                    key = line.pop(0)
                    if key==')':
                        # end of 'actions' block
                        return None
                    if is_array:
                        assert(len(line)==3)
                        value = np.array([float(val) for val in line])
                    else:
                        assert(len(line)==1)
                        try:
                            value = float(line[0])
                        except ValueError:
                            value = line[0]
                    inputs[key] = value
            return inputs

        # read all blocks
        N = 0
        while True:
            inputs = read_block()
            if inputs is None:
                break
            N += 1
            if verbose:
                print('{:d}: {} source {}'.format(
                        N,inputs['action'],inputs['source']))
            if (IDs is not None) and (N not in IDs):
                if verbose: print('  skipped')
                continue
            plot_info = plot_kwargs.copy()
            if N==1:
                name = '{:s} : {:s}'.format(
                        os.path.split(fpath)[-1], inputs['source'])
                plot_info['label'] = name
            funcname = 'plot_' + inputs['source']
            try:
                plotfun = getattr(sys.modules[__name__],funcname)
            except AttributeError:
                print(funcname,'not available')
            else:
                plotfun(plane,plot_kwargs=plot_info,**inputs)

def _plot_poly(*args,**kwargs):
    args = list(args)
    args.append(args[0])
    pts = np.array(args)
    plt.plot(pts[:,0],pts[:,1],**kwargs)
    plt.axis('equal')

def plot_rotatedBoxToCell(plane,plot_kwargs={},**kwargs):
    origin = np.array(kwargs['origin'])
    if plane=='xy':
        # horizontal slice
        origin = origin[[0,1]]
        r1 = np.array(kwargs['i'])[[0,1]]
        r2 = np.array(kwargs['j'])[[0,1]]
    elif plane=='xz':
        # vertical slice
        origin += np.array(kwargs['j'])/2 # slice through centroid of box
        origin = origin[[0,2]]
        r1 = np.array(kwargs['i'])[[0,2]]
        r2 = np.array(kwargs['k'])[[0,2]]
    else:
        print('unknown plane orientation:',plane)
    LL = origin
    LR = origin + r1
    UR = origin + r1 + r2
    UL = origin + r2
    _plot_poly(LL,LR,UR,UL,**plot_kwargs)

def plot_cylinderToCell(plane,plot_kwargs={},**kwargs):
    p1 = np.array(kwargs['p1'])
    p2 = np.array(kwargs['p2'])
    R = kwargs['radius']
    if plane=='xy':
        # horizontal slice
        p1 = p1[[0,1]]
        p2 = p2[[0,1]]
    elif plane=='xz':
        # vertical slice
        p1 = p1[[0,2]]
        p2 = p2[[0,2]]
    else:
        print('unknown plane orientation:',plane)
    dp = p2 - p1 
    ang = np.arctan2(dp[1],dp[0])
    dr = np.array([-dp[1],dp[0]])
    dr /= np.sqrt(dr.dot(dr))
    LL = p1 - dr*R
    LR = p2 - dr*R
    UR = p2 + dr*R
    UL = p1 + dr*R
    #plt.plot([p1[0],p2[0]], [p1[1],p2[1]], 'k--')
    _plot_poly(LL,LR,UR,UL,**plot_kwargs)


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
    """Class for generating topoSetDict files. To visualize existing 
    topoSetDicts, use the plot_toposet function.
    """
    source_types = ['box','cylinder']
    defaults = dict(
        rotation=0.0,
        upstream=5.0, downstream=10.0,
        cyl_upstream=0.5, cyl_downstream=0.5,
        width=3.0, height=3.0,
        xbuffer_upstream=1.0,
        xbuffer_downstream=1.0,
        ybuffer=1.0, zbuffer=1.0,
        rbuffer=0.5,
        zoffset=0.0
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
        self.cyl_upstream = []
        self.cyl_downstream = []
        self.width = []
        self.height = []
        self.zoffset = []
        self.xbuffer_upstream = []
        self.xbuffer_downstream = []
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
        upstream, downstream : float
            Distance (in diameters) upstream/downstream of the turbine
            where the innermost refinement box starts.
        cyl_upstream, cyl_downstream : float
            Distance (in diameters) upstream/downstream of the turbine
            where the innermost refinement cylinder starts.
        width : float
            The overall width (in diameters) of the inner refinement
            box; need to account for horizontal wake meandering
            downstream.
        height : float
            The overall height (in diameters) of the inner refinement
            box; need to account for vertical wake meandering
            downstream.
        xbuffer : float or list-like
        xbuffer_upstream : float or list-like
        xbuffer_downstream : float or list-like
            Size of buffer region (in diameters) between refinement
            levels in the upstream/downstream directions; set 'xbuffer'
            to use the same value for upstream and downstream or set
            'xbuffer_upstream' and 'xbuffer_downstream' separately--for
            box and cylinder refinement
        ybuffer : float or list-like
            Size of buffer region (in diameters) between refinement
            levels in the lateral directions--for box and cylinder 
            refinement
        zbuffer : float or list-like
            Size of buffer region (in diameters) between refinement
            levels in the vertical directions--for box refinement only
        rbuffer : float
            Used to set size of cylindrical refinement region; cylinder
            level i has diameter (1 + (i+1)*rbuffer)*D for
            i=0,1,...
        zoffset : float
            For box refinement region, the height/elevation of the lower
            surface offset from the turbine base location z value; set
            to < 0 to capture regions on windward/leeward regions if a
            turbine is situated on a hill. [m]
        """
        # special treatment of xbuffer inputs
        if 'xbuffer' in kwargs:
            if not 'xbuffer_upstream' in kwargs:
                kwargs['xbuffer_upstream'] = kwargs['xbuffer']
            if not 'xbuffer_downstream' in kwargs:
                kwargs['xbuffer_downstream'] = kwargs['xbuffer']
        # now plug in all the remaining values
        for key,defval in self.defaults.items():
            val = kwargs.get(key, defval)
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
        refinement sources. Turbine-specific refinement parameters may
        be specified as keyword arguments, otherwise the default values
        set with setup() are used.

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
        self.base_location.append(np.array(base_location,dtype=float))
        self.diameter.append(D)
        if zhub is None:
            zhub = D
        self.zhub.append(zhub)
        def get_param(param):
            return kwargs.get(param, self.refinement[param])
        self.rotation.append(get_param('rotation'))
        self.upstream.append(get_param('upstream'))
        self.downstream.append(get_param('downstream'))
        self.cyl_upstream.append(get_param('cyl_upstream'))
        self.cyl_downstream.append(get_param('cyl_downstream'))
        self.width.append(get_param('width'))
        self.height.append(get_param('height'))
        self.zoffset.append(get_param('zoffset'))
        self.xbuffer_upstream.append(get_param('xbuffer_upstream'))
        self.xbuffer_downstream.append(get_param('xbuffer_downstream'))
        self.ybuffer.append(get_param('ybuffer'))
        self.zbuffer.append(get_param('zbuffer'))
        self.rbuffer.append(get_param('rbuffer'))


    def estimate_mesh_size(self,initial_size=None,ds0=10.0):
        """Estimate mesh size from the initial mesh size (e.g., from
        blockMesh) and an initial uniform mesh size.

        NOTE: this is experimental and not very accurate!
        """
        if initial_size is None:
            raise ValueError('specify initial cell count or list of dimensions')
        if hasattr(initial_size,'__iter__'):
            initial_size = int(np.prod(initial_size))
            print('calculated initial cell count: {:d}'.format(initial_size))
        Nlevels = len(self.sources)
        vol = np.zeros(Nlevels)
        for ilevel in range(Nlevels-1,-1,-1):
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
                    xbuff_u = self._get_refinement_buffer(self.xbuffer_upstream[iturb],efflevel) * Lref
                    xbuff_d = self._get_refinement_buffer(self.xbuffer_downstream[iturb],efflevel) * Lref
                    ybuff = self._get_refinement_buffer(self.ybuffer[iturb],efflevel) * Lref
                    zbuff = self._get_refinement_buffer(self.zbuffer[iturb],efflevel) * Lref
                    length += xbuff_u + xbuff_d
                    width += 2*ybuff
                    height += zbuff
                    vol[ilevel] += length*width*height
            elif sourcename == 'cylinder':
                print('{:d}: turbine cylinder regions {:d}'.format(ilevel,efflevel))
                for iturb in range(len(self.base_location)):
                    Lref = self.diameter[iturb]
                    R = (1.0 + (efflevel+1)*self.rbuffer[iturb]) * Lref/2
                    upstream = self.upstream[iturb] * Lref
                    downstream = self.downstream[iturb] * Lref
                    xbuff_u = self._get_refinement_buffer(self.xbuffer_upstream[iturb],efflevel) * Lref
                    xbuff_d = self._get_refinement_buffer(self.xbuffer_downstream[iturb],efflevel) * Lref
                    length = upstream + downstream + xbuff_u + xbuff_d
                    vol[ilevel] += length * np.pi*R**2
        lastcount = initial_size
        ds = float(ds0)
        for ilevel in range(Nlevels-1,-1,-1):
            print('volume, ds : {:g} {:f}'.format(vol[ilevel],ds))
            approx_cells = int(vol[ilevel] / ds**3)
            print('approx number of selected cells : {:d}'.format(approx_cells))
            newcount = lastcount + 7*approx_cells
            print('refinement level {:d} : increase from {:d} to {:d}'.format(
                    ilevel+1, lastcount, newcount))
            lastcount = newcount
            ds /= 2


    def plot(self,plane='xy',turbines=None):
        """Visualize locations of turbines. If 'turbines' is None, all
        are plotted, otherwise a list of one-indexed IDs should be
        specified.

        TODO: call plot_* functions here too
        """
        R = np.array(self.diameter) / 2
        if plane=='xy':
            locations = np.array([ loc[[0,1]] for loc in self.base_location ])
            ang = np.deg2rad(self.rotation)
            xoff = R * np.sin(ang)
            yoff = R * np.cos(ang)
            xrotor = [ [xval+xoff,xval-xoff] for xval in locations[:,0] ]
            yrotor = [ [yval-yoff,yval+yoff]
                       for iturb,yval in enumerate(locations[:,1]) ]
        elif plane=='xz':
            locations = np.array([ loc[[0,2]] for loc in self.base_location ])
            locations[:,1] += self.zhub
            xrotor = [ [xval,xval] for xval in locations[:,0] ]
            yrotor = [ [zval-R[iturb],zval+R[iturb]]
                       for iturb,zval in enumerate(locations[:,1]) ]
        else:
            print('unknown plane orientation:',plane)
        if turbines is None:
            turbines = np.arange(len(locations))
        else:
            turbines = [iturb-1 for iturb in turbines]
        for iturb in turbines:
            loc = locations[iturb]
            plt.plot(loc[0],loc[1], 'ko', markersize=5, markerfacecolor='k')
            plt.plot(xrotor[iturb],yrotor[iturb], 'k-')
            plt.text(loc[0],loc[1], '{:d}'.format(iturb+1), fontsize='xx-large')


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
        # bkg_ybuffer, and bkg_zbuffer which are all dimensional
        # parameters.
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
        ang = np.deg2rad(self.bkg_rotation[ibkg])
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


    def _get_refinement_buffer(self,buff,ilevel):
        """*buffer parameters may be set to a scalar or a list-like
        object. This returns the size of the buffer region up to the 
        specified refinement nest level.
        """
        if hasattr(buff,'__iter__'):
            return np.sum(buff[:ilevel])
        else:
            return ilevel*buff


    def _write_box(self,iturb,ilevel):
        """Write out topoSetDict for a turbine refinement box.

        Depends on D, upstream, downstream, width, height,
        xbuffer_upstream, xbuffer_downstream, ybuffer, and zbuffer.
        """
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
        zoffset = self.zoffset[iturb]
        #xbuff = ilevel * self.xbuffer[iturb] * Lref
        xbuff_u = self._get_refinement_buffer(self.xbuffer_upstream[iturb],ilevel) * Lref
        xbuff_d = self._get_refinement_buffer(self.xbuffer_downstream[iturb],ilevel) * Lref
        ybuff = self._get_refinement_buffer(self.ybuffer[iturb],ilevel) * Lref
        zbuff = self._get_refinement_buffer(self.zbuffer[iturb],ilevel) * Lref
        base = self.base_location[iturb]
        ang = np.deg2rad(self.rotation[iturb])
        # origin (x,y,z)
        x0 = -upstream - xbuff_u
        y0 = -0.5*width - ybuff
        x = x0 * np.cos(ang) - y0 * np.sin(ang) + base[0]
        y = x0 * np.sin(ang) + y0 * np.cos(ang) + base[1]
        z = base[2] + zoffset
        # box dimensions
        #ix = (length + 2*xbuff) * np.cos(ang)
        #iy = (length + 2*xbuff) * np.sin(ang)
        ix = (length + xbuff_u + xbuff_d) * np.cos(ang)
        iy = (length + xbuff_u + xbuff_d) * np.sin(ang)
        iz = 0.0
        jx = -(width + 2*ybuff) * np.sin(ang)
        jy =  (width + 2*ybuff) * np.cos(ang)
        jz = 0.0
        kx = 0.0
        ky = 0.0
        kz = height - zoffset + zbuff
        return template.format(action=action,
                x0=x, y0=y, z0=z,
                ix=ix, iy=iy, iz=iz,
                jx=jx, jy=jy, jz=jz,
                kx=kx, ky=ky, kz=kz)


    def _write_cylinder(self,iturb,ilevel):
        """Write out topoSetDict for a turbine refinement cylinder.

        Depends on D, radial_buffer, cyl_upstream, cyl_downstream,
        xbuffer_upstream, and xbuffer_downstream.
        """
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
        upstream = self.cyl_upstream[iturb] * Lref
        downstream = self.cyl_downstream[iturb] * Lref
        #xbuff = ilevel * self.xbuffer[iturb] * Lref
        xbuff_u = self._get_refinement_buffer(self.xbuffer_upstream[iturb],ilevel) * Lref
        xbuff_d = self._get_refinement_buffer(self.xbuffer_downstream[iturb],ilevel) * Lref
        base = self.base_location[iturb]
        ang = np.deg2rad(self.rotation[iturb])
        # upstream point
        xu = -upstream - xbuff_u
        x1 = xu * np.cos(ang) + base[0]
        y1 = xu * np.sin(ang) + base[1]
        xd = downstream + xbuff_d
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

