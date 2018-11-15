#
# Utilities for setting up OpenFOAM sampling methods
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
# Note:
# - function object ref: OpenFOAM-2.4.x/src/postProcessing/functionObjects/IO/controlDict
#
from __future__ import print_function
import sys
import numpy as np

class SampleSet(list):
    """A container for "set" sampling objects"""

    header = """{name:s}
{{
    type                sets;
    functionObjectLibs ("libsampling.so" "libuserfileFormats.so");
    enabled             true;
    interpolationScheme {interpolationScheme:s};
    outputControl       {outputControl:s};
    {intervalName:s}      {outputInterval:g};
    setFormat           {setFormat:s};
    fields
    (
        {fields:s}
    );

    sets
    (
"""
    footer = """    );
}
"""

    def __init__(self,name,*args,**kwargs):
        self.name = name
        self.fields = kwargs.pop('fields',['U'])
        self.interpolationScheme = kwargs.pop('interpolationScheme','cellPoint')
        self.outputControl = kwargs.pop('outputControl','timeStep')
        self.outputInterval = kwargs.pop('outputInterval',1)
        self.setFormat = kwargs.pop('setFormat','NEED_TO_SPECIFY')
        if sys.version_info < (3, 0):
            super(SampleSet, self).__init__(*args,**kwargs)
        else:
            super().__init__(*args,**kwargs)

    def __repr__(self):
        N = len(self)
        if N == 0:
            s = 'empty SampleSet'
        else:
            s = '{:d} sets in sampleSet\n'.format(N)
        for i,sampleset in enumerate(self):
            s += '{:d}: {:s}\n'.format(i,str(sampleset))
        return s

    def write(self,fpath):
        """Write sampling object definition"""
        if self.outputControl in ['timeStep','outputTime']:
            intervalName = 'outputInterval'
            intervalValue = int(self.outputInterval)
        else:
            intervalName = 'writeInterval '
            intervalValue = self.outputInterval
        fields = '\n        '.join(self.fields)
        with open(fpath,'w') as f:
            f.write(
                    self.header.format(
                        name=self.name,
                        interpolationScheme=self.interpolationScheme,
                        outputControl=self.outputControl,
                        outputInterval=intervalValue,
                        intervalName=intervalName,
                        setFormat=self.setFormat,
                        fields=fields
                        )
                   )
            for sampleset in self:
                sampleset.write(f)
            f.write(self.footer)

    def plotxy(self,ax=None,label=True,**kwargs):
        """Plot all bounding boxes from the sampleSet in an xy-plane"""
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        if ax is None:
            fig,ax = plt.subplots()
        for i,sampleset in enumerate(self):
            LL = sampleset.bbox0[:2] # lower-left
            UR = sampleset.bbox1[:2] # upper-right
            center = (LL+UR) / 2
            extent = UR - LL
            rect = patches.Rectangle(LL,width=extent[0],height=extent[1],
                                     fill=False, **kwargs)
            ax.add_patch(rect)
            if label:
                if not (label==True):
                    label = label.format(i+1)
                ax.text(center[0],center[1],sampleset.name,
                        horizontalalignment='center',
                        verticalalignment='center',
                        **kwargs)
        

class Set(object):
    """Base class for sampling sets"""
    def __init__(self,settype,template,perturb=0.001):
        self.name = 'NEED_TO_SPECIFY'
        self.type = settype
        self.template = template
        self.perturb = perturb
        self.params = {'type': settype}

    def write(self,f):
        """Write a single sampling set definition"""
        params = dict(name=self.name,type=self.type)
        for key,val in self.params.items():
            params[key] = val
        f.write(self.template.format(**params))


class Array(Set):
    """OpenFOAM "array" style sampling, e.g., for FAST.Farm inflow

    Reference: https://cpp.openfoam.org/v5/arraySet_8C_source.html
    """

    template = """        {name:s}
        {{
            type                {type:s};
            axis                {axis:s};
            origin              ({x0:.3f} {y0:.3f} {z0:.3f});
            coordinateRotation
            {{
                type            {rotation:s};
                e1              {e1:s};
                e2              {e2:s};
            }}
            pointsDensity       ({pointsDensity[0]:d} {pointsDensity[1]:d} {pointsDensity[2]:d});
            spanBox             ({spanBox[0]:.3f} {spanBox[1]:.3f} {spanBox[2]:.3f});
        }}
"""

    def __init__(self,name,origin,
                      pointsDensity=None, spanBox=None,
                      dimensions=None, spacing=None,
                      axis='xyz',
                      rotation='axesRotation',e1=(1,0,0),e2=(0,1,0),
                      **kwargs):
        """Generate sampling input file for OpenFOAM given the 'origin' and two
        of the following:
        - 'pointsDensity', 'spanBox' (direct inputs to OpenFOAM)
        - 'dimensions', 'spacing' (description of sampling grid)
        """
        if sys.version_info < (3, 0):
            super(self.__class__, self).__init__('array',self.template,**kwargs)
        else:
            super().__init__('array',self.template,**kwargs)

        assert( (pointsDensity and spanBox) or (dimensions and spacing) )
        self.name = name
        origin = np.array(origin)

        # calculate dimensions/spacing or pointsDensity/spanBox
        if (pointsDensity is not None) and (spanBox is not None):
            self.params['pointsDensity'] = np.array(pointsDensity, dtype=int)
            self.params['spanBox'] = np.array(spanBox)
            self.spacing = self.params['spanBox'] / (self.params['pointsDensity'] + 1)
        else:
            # calculate origin and spanBox
            self.spacing = np.array(spacing)
            self.params['pointsDensity'] = np.array(dimensions, dtype=int)
            self.params['spanBox'] = self.spacing * (self.params['pointsDensity'] + 1)
            origin -= self.spacing

        # note: per the OpenFOAM sampling definition, the specified origin
        # is _not_ actually a sampling point
        self.bbox0 = origin + self.spacing
        self.bbox1 = origin + self.params['pointsDensity']*self.spacing

        # perturb sampling points to ensure we don't end up on a face or edge
        # as this sometimes causes the sampling to fail (an OpenFOAM quirk!)
        self.params['x0'] = origin[0] + self.perturb
        self.params['y0'] = origin[1] + self.perturb
        self.params['z0'] = origin[2] + self.perturb

        # optional parameters
        self.params['axis'] = axis
        self.params['rotation'] = rotation
        self.params['e1'] = str(e1).replace(',','')
        self.params['e2'] = str(e2).replace(',','')

    def __repr__(self):
        s = '"{:s}" array\n'.format(self.name)
        s += '  bounding box ({:s}, {:s})\n'.format(str(self.bbox0),
                                                     str(self.bbox1))
        s += '  dimensions {:s}\n'.format(str(self.params['pointsDensity']))
        s += '  spacing {:s}'.format(str(self.spacing))
        return s


class Probes(list):
    """General point sampling"""

    header = """{name:s}
{{
    type                probes;
    functionObjectLibs ("libsampling.so");
    name                {name:s};
    outputControl       {outputControl:s};
    {intervalName:s}      {outputInterval:g};
    interpolationScheme {interpolationScheme:s};
    fields
    (
        {fields:s}
    );

    probeLocations
    (
"""
    footer = """    );
}
"""

    def __init__(self,name,*args,**kwargs):
        """Probe sampling dictionary, to be included in controlDict.functions

        Locations may be specified as a list of coordinates or as a 2D array.

        Keyword arguments
        -----------------
        fields : list, optional
            Fields to sample (default: U, T)
        outputControl : str, optional
            'timeStep', 'runTime', or 'adjustableRunTime' (default: timeStep)
        outputInterval : int or float, optional
            Number of steps or sampling period, for timestep and runtime
            sampling, respectively (default: 1)
        interpolationScheme : str, optional
            'cell', 'cellPoint' linear (default: cell)
        perturb : float, optional
            Shift sampling location by a small amount to help prevent the probe
            from landing on an edge or face; under some circumstances, OpenFOAM
            sampling may fail if this is the case.
        """
        self.name = name
        self.fields = kwargs.pop('fields',['U','T'])
        self.outputControl = kwargs.pop('outputControl','timeStep')
        self.outputInterval = kwargs.pop('outputInterval',1)
        self.interpolationScheme = kwargs.pop('interpolationScheme','cell')
        self.perturb = kwargs.pop('perturb',0.001)
        if sys.version_info < (3, 0):
            super(Probes, self).__init__(*args,**kwargs)
        else:
            super().__init__(*args,**kwargs)

    def __repr__(self):
        return '{:s} ({:d} probes)'.format(self.name, len(self))

    def write(self,fpath=None):
        """Write a point sampling object definition"""
        if fpath is None:
            fpath = self.name
        fields = '\n        '.join(self.fields)
        if self.outputControl in ['timeStep','outputTime']:
            intervalName = 'outputInterval'
            intervalValue = int(self.outputInterval)
        else:
            intervalName = 'writeInterval '
            intervalValue = self.outputInterval
        with open(fpath,'w') as f:
            f.write(
                    self.header.format(
                        name=self.name,
                        outputControl=self.outputControl,
                        outputInterval=intervalValue,
                        intervalName=intervalName,
                        interpolationScheme=self.interpolationScheme,
                        fields=fields
                        )
                   )
            locations = np.array(self)
            for i in range(len(self)):
                loc = locations[i,:]
                f.write('        ({:f} {:f} {:f})\n'.format(loc[0]+self.perturb,
                                                            loc[1]+self.perturb,
                                                            loc[2]+self.perturb))
            f.write(self.footer)
        print('Wrote '+fpath)


class VirtualMetMast(Probes):
    def __init__(self,name,base,*args,**kwargs):
        """Create a virtual met mast described by a Probes object

        base: list-like
            xyz coordinates of the base of the met mast
        heights: list-like, optional
            If specified, calls add_heights
        """
        heights = kwargs.pop('heights',None)
        if sys.version_info < (3, 0):
            super(VirtualMetMast, self).__init__(name,*args,**kwargs)
        else:
            super().__init__(name,*args,**kwargs)
        self.base = np.array(base)
        if heights is not None:
            self.add_heights(heights)

    def add_heights(self,z):
        for zi in z:
            self.append([self.base[0],self.base[1],self.base[2]+zi])
