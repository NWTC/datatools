header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.0                                   |
|   \\\\  /    A nd           | Web:      http://www.OpenFOAM.org               |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbineArrayProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
globalProperties
{{
    outputControl       {outputControl:s};
    outputInterval      {outputInterval:g};
}}

"""

turbine_definition = """{name:s}
{{
    turbineType                      "{turbineType:s}";
    includeNacelle                    {includeNacelle:s};
    includeTower                      {includeTower:s};
    baseLocation                     ({baseLocation[0]:g} {baseLocation[1]:g} {baseLocation[2]:g});
    numBladePoints                    {numBladePoints:d};
    numNacellePoints                  {numNacellePoints:d};
    numTowerPoints                    {numTowerPoints:d};
    bladePointDistType               "{bladePointDistType:s}";
    nacellePointDistType             "{nacellePointDistType:s}";
    towerPointDistType               "{towerPointDistType:s}";
    bladeSearchCellMethod            "{bladeSearchCellMethod:s}";
    bladeActuatorPointInterpType     "integral";
    nacelleActuatorPointInterpType   "linear";
    towerActuatorPointInterpType     "linear";
    actuatorUpdateType               "{actuatorUpdateType:s}";
    velocityDragCorrType             "{velocityDragCorrType:s}";
    bladeForceProjectionType         "{bladeForceProjectionType:s}";
    nacelleForceProjectionType       "{nacelleForceProjectionType:s}";
    towerForceProjectionType         "{towerForceProjectionType:s}";
    bladeForceProjectionDirection    "{bladeForceProjectionDirection:s}";
    bladeEpsilon                     ({bladeEpsilon[0]:f} {bladeEpsilon[1]:f} {bladeEpsilon[2]:f});
    nacelleEpsilon                   ({nacelleEpsilon[0]:f} {nacelleEpsilon[1]:f} {nacelleEpsilon[2]:f});
    towerEpsilon                     ({towerEpsilon[0]:f} {towerEpsilon[1]:f} {towerEpsilon[2]:f});
    nacelleSampleDistance             {nacelleSampleDistance:f};
    towerSampleDistance               {towerSampleDistance:f};
    tipRootLossCorrType              "{tipRootLossCorrType:s}";
    rotationDir                      "{rotationDir:s}";
    Azimuth                           {azimuth:f};
    RotSpeed                          {rotorSpeed:f};
    TorqueGen                         {torqueGen:f};
    Pitch                             {pitch:f};
    NacYaw                            {nacelleYaw:f};
    fluidDensity                      {fluidDensity:f};
}}

"""

defaults = dict(
    includeNacelle='false',
    includeTower='false',
    numNacellePoints=10,
    numTowerPoints=80,
    bladePointDistType='uniform',  # (uniform | ...)
    nacellePointDistType='uniform',  # (uniform | ...)
    towerPointDistType='uniform',  # (uniform | ...)
    bladeSearchCellMethod='disk',  # (sphere | disk)
    bladeActuatorPointInterpType='integral',  # (integral | linear | cellCenter)
    nacelleActuatorPointInterpType='linear',  # (linear | cellCenter)
    towerActuatorPointInterpType='linear',  # (linear | cellCenter)
    actuatorUpdateType='oldPosition',  # (oldPosition | newPosition)
    velocityDragCorrType='none',  # (none | Martinez)
    # bladeForceProjectionType: (uniformGaussian | generalizedGaussian | generalizedGaussian2D | variableGaussianUserDef | variableUniformGaussianChord | chordThicknessGaussian | chordThicknessGaussian2D)
    bladeForceProjectionType='uniformGaussian',
    # nacelleForceProjectionType: (uniformGaussian | diskGaussian | advanced*)
    nacelleForceProjectionType='diskGaussian',
    # towerForceProjectionType: (uniformGaussian | diskGaussian | ringGaussian | advanced*)
    towerForceProjectionType='diskGaussian',
    # bladeForceProjectionDirection: (localVelocityAligned | localVelocityAlignedCorrected | sampleVelocityAligned)
    bladeForceProjectionDirection='localVelocityAligned',
    nacelleEpsilon=(0,0,0),
    towerEpsilon=(0,0,0),
    nacelleSampleDistance=1.0,
    towerSampleDistance=1.0,
    tipRootLossCorrType='Glauert',  # (none | Glauert)
    rotationDir='cw',  # (cw | ccw)
    fluidDensity=1.2,
)


class TurbineArrayProperties(object):
    """Class for creating constant/turbineArrayProperties input files"""
    def __init__(self,outputControl='timeStep',outputInterval=1,
                 **kwargs):
        """Optional turbine properties are specified as keyword
        arguments. See turbineArrayProperties.defaults for default
        parameters.
        """
        self.outputControl = outputControl;
        self.outputInterval = outputInterval;
        self.properties = defaults.copy()
        self.turbines = []


    def add_turbine(self,turbineType,
                    baseLocation,numBladePoints,bladeEpsilon,
                    name=None,
                    azimuth=0.0,
                    rotorSpeed=0.0,
                    torqueGen=0.0,
                    pitch=0.0,
                    nacelleYaw=270.0,
                    **kwargs):
        """Add turbine at 'baseLocation' with the specified number of
        actuator points and spreading parameters (epsilon). 
        'turbineType' should correspond to a turbine definition in
        constant/turbineProperties/<turbineType>.

        Operating parameters
        --------------------
        name : str, optional
            if None, then defaults 'turbine${i}' where i = 0..Nturb-1
        azimuth : float, optional
            Blade 1 azimuth [deg]
        rotorSpeed : float
            Rotor RPM
        torqueGen : float, optional
            commanded generator torque, used to calculate speed if
            generator torque or blade pitch control is turned on [N-m]
        pitch : float
            blade pitch [deg]
        nacelleYaw: float
            Nacelle yaw angle, as a compass direction [deg]

        Default turbine properties may be overridden using keyword
        arguments.
        """
        d = self.properties.copy()
        if name  is None:
            name = 'turbine{:d}'.format(len(self.turbines))
        d['name'] = name
        d['turbineType'] = turbineType
        d['baseLocation'] = baseLocation
        d['numBladePoints'] = numBladePoints
        d['bladeEpsilon'] = bladeEpsilon
        d['azimuth'] = azimuth
        d['rotorSpeed'] = rotorSpeed
        d['torqueGen'] = torqueGen
        d['pitch'] = pitch
        d['nacelleYaw'] = nacelleYaw
        for key,val in kwargs.items():
            if not key in d:
                print("Warning: specified parameter '{:s}' not found")
            elif isinstance(val,bool):
                d[key] = str(val).lower()
            else:
                d[key] = val
        self.turbines.append(d)


    def write(self,fpath='turbineArrayProperties'):
        """Write OpenFOAM dictionary constant/turbineArrayProperties"""
        with open(fpath,'w') as f:
            f.write(header.format(outputControl=self.outputControl,
                                  outputInterval=self.outputInterval))
            for iturb, d in enumerate(self.turbines):
                f.write(turbine_definition.format(**d))
        print('Wrote',fpath)

