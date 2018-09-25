#!/usr/bin/env python
#
# Script to convert all averaging data to a pandas dataframe
# written by Eliot Quon (eliot.quon@nrel.gov)
#
import sys
from datatools.SOWFA.postProcessing.averaging import PlanarAverages
avg = PlanarAverages( *sys.argv[1:] )
avg.to_csv('averaging.csv',fields='all')
