#
# Utilities for setting up OpenFOAM sampling methods
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import sys
import numpy as np

class SampleSet(list):
    """A container for "set" sampling objects"""

    header = """{name:s}
{{
    type                sets;
    functionObjectLibs ("libsampling.so");
    enabled             true;
    interpolationScheme {interpolationScheme:s};
    outputControl       timeStep;
    outputInterval      {outputInterval:d};
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
        fields = '\n        '.join(self.fields)
        with open(fpath,'w') as f:
            f.write(
                    self.header.format(
                        name=self.name,
                        interpolationScheme=self.interpolationScheme,
                        outputInterval=self.outputInterval,
                        setFormat=self.setFormat,
                        fields=fields
                        )
                   )
            for sampleset in self:
                sampleset.write(f)
            f.write(self.footer)
        

class Set(object):
    """Base class for sampling sets"""
    def __init__(self,settype,template,perturb=0.001):
        self.name = 'NEED_TO_SPECIFY'
        self.type = settype
        self.template = template
        self.perturb = perturb
        self.params = {'type': settype}

    def write(self,f):
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
            pointsDensity       (12 12 13);
            spanBox             (130 130 140);
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

