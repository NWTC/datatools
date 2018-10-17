#!/usr/bin/env python
import numpy as np
import datatools.SOWFA.constant.boundaryData as bd

# TODO: setup 1D times array

# TODO: setup 1D x,y,z arrays
patch = bd.CartesianPatch(x,y,z,dpath='boundaryData',name='west')
patch.write_points()

# TODO: generate U,T arrays with shape==(time, height[, direction])
patch.write_profiles(times,z,U=U,T=T)

