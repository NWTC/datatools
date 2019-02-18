#!/usr/bin/env python
#
# For creating symlinks to output from the OpenFOAM sample utility
#
# Expected directory structure:
#   postProcessing/surfaces/<time>/<var>_<sampleName>.vtk
# Result:
#   postProcessing/surfaces/<sampleName>/<var>_<time>.vtk
#
# If the -nameFirst argument is passed, then the expected output files
# should be named <sampleName>_<var>.vtk
#

from __future__ import print_function
import os
import sys
import argparse

# workarounds:
extMapping = dict(xy='xyz')
expectedVars = ['U','T','k','p','p_rgh']

# parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('srcdir',nargs='?',default=os.getcwd(),help='alternate source directory')
parser.add_argument('--vars',nargs='+',default=expectedVars,
                    help='specify expected variable names')
parser.add_argument('-v','--verbose',action='store_true',
                    help='verbose output (for debugging)')
args = parser.parse_args()

verbose = args.verbose

# find time subdirectories
dirlist = []
timesteps = []
srcdir = args.srcdir
dirs = os.walk(srcdir).next()[1]
for d in dirs:
    try: 
        step = float(d) # need this to verify this is a time-step dir!
    except ValueError:
        pass
    else:
        dirlist.append(os.path.join(srcdir,d))
        timesteps.append(step)
if len(dirlist) == 0:
    sys.exit('No time subdirectories found in '+str(dirs))

# set up search strings for the expected variables with underscore vars first
varlist = args.vars
expectedVars = [v for v in varlist if v.__contains__('_')] \
             + [v for v in varlist if not v.__contains__('_')]
if verbose:
    print('Expected variables:',expectedVars)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

# now, find sample and variable name(s)
sampleNames = []
fullVarNames = [] # consisting of one or more variable names separated by '_'
varNames = []
extNames = []
ext = None
nameFirst = None
for timestep_dir in dirlist:
    if verbose:
        print('Processing', timestep_dir)
    filelist = [ f for f in os.listdir(timestep_dir)
                 if os.path.isfile(os.path.join(timestep_dir,f)) ]
    for f in filelist:
        if f.startswith('.'):
            continue
        name,ext = os.path.splitext(f)
        basename = name
        ext = ext[1:] # drop leading '.' from extension
        newvars = []
        for var in expectedVars:
            if var in name:
                newvars.append(var)
                name = name.replace(var,'')
        if name.endswith('_'):
            nameFirst = True
            name = name.rstrip('_')
        elif name.startswith('_'):
            nameFirst = False
            name = name.lstrip('_')
        else:
            name = name.replace('_','')
        if verbose:
            print('  {:s}\t(name={:s}, var={:s}, ext={:s})'.format(f,name,','.join(newvars),ext))
        if name=='':
            name = 'timeSeries'
        else:
            if nameFirst:
                fullvarname = basename.replace(name+'_','')
            else:
                fullvarname = basename.replace('_'+name,'')
        if not name in sampleNames:
            sampleNames.append(name)
            if not os.path.exists(name): os.makedirs(name)
        for var in newvars:
            if var not in varNames:
                varNames.append(var)
                fullVarNames.append(fullvarname)
        if not ext in extNames:
            extNames.append(ext)
if ext is None:
    sys.exit('No vtk files found in '+str(dirlist))
if not len(extNames)==1:
    sys.exit('Don''t know how to handle different extensions',extNames)
if ext in extMapping:
    extNew = extMapping[ext]
else:
    extNew = ext

# create symlinks
indices = sorted(range(len(timesteps)), key=lambda k: timesteps[k])
for sample in sampleNames:
    for var,fullvar in zip(varNames,fullVarNames):
        for i in range(len(timesteps)):
            idx = indices[i]
            dname = dirlist[idx]#.split()
            if sample=='timeSeries':
                src = os.path.join( os.getcwd(), dname, var+'.'+ext )
            else:
                if nameFirst:
                    src = os.path.join( os.getcwd(), dname, sample+'_'+fullvar+'.'+ext )
                else:
                    src = os.path.join( os.getcwd(), dname, fullvar+'_'+sample+'.'+ext )
            dest = os.path.join(sample, '{:s}_{:d}.{:s}'.format(var,i,extNew))
            if verbose:
                print(dest,'-->',src)
            try:
                os.symlink(src,dest)
            except OSError:
                pass

print('sample names: ',sampleNames)
print('field names: ',varNames)
print('full field names: ',fullVarNames)
print('time steps: ',len(timesteps),'[',timesteps[indices[0]],'...',timesteps[indices[-1]],']')

