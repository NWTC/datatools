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

parser = argparse.ArgumentParser()
parser.add_argument('srcdir',nargs='?',default=os.getcwd(),help='alternate source directory')
parser.add_argument('--nameFirst',action='store_false',
                    help='parse output files assuming the sampleName preceeds the variable name(s)')
parser.add_argument('-v','--verbose',action='store_false',
                    help='verbose output (for debugging)')
args = parser.parse_args()

#verbose = False
verbose = True # for debug

dirlist = []
timesteps = []

if len(sys.argv) > 1:
    srcdir = sys.argv[1]
else:
    srcdir = '.'

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

extMapping = dict(xy='xyz')

underscoreNames = ['p_rgh']

def tname(tval):
    # note: paraview doesn't seem to handle floats well...
    #return '%d' % (tval*10)
    #return '%d' % (tval*1000)
    return '%d' % (tval)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

sampleNames = []
varNames = []
extNames = []
ext = None
for timestep_dir in dirlist:
    if verbose:
        print('Processing', timestep_dir)
    filelist = [ f for f in os.listdir(timestep_dir)
                 if os.path.isfile(os.path.join(timestep_dir,f)) ]
    for f in filelist:
        if f.startswith('.'):
            continue
        basename,ext = os.path.splitext(f)
        ext = ext[1:] # drop leading '.' from extension
        origname = None
        for uname in underscoreNames:
            if uname in basename:
                origname = uname
                basename = basename.replace(uname,'TEMP')
        basesplit = basename.split('_')
        var = basesplit[0]
        if origname is not None:
            var = var.replace('TEMP',origname)
        name = '_'.join(basesplit[1:])
        if verbose:
            print('  {:s}\t(name={:s}, var={:s}, ext={:s})'.format(f,name,var,ext))
        if name=='':
            name = 'timeSeries'
        if not name in sampleNames:
            sampleNames.append(name)
            if not os.path.exists(name): os.makedirs(name)
        if not var in varNames:
            varNames.append(var)
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

indices = sorted(range(len(timesteps)), key=lambda k: timesteps[k])
for sample in sampleNames:
    for var in varNames:
        for i in range(len(timesteps)):
            idx = indices[i]
            dname = dirlist[idx]#.split()
            if sample=='timeSeries':
                src = os.path.join( os.getcwd(), dname, var+'.'+ext )
                dest = sample + os.sep + '%s_%s.%s' % (var,i,extNew)
            else:
                src = os.path.join( os.getcwd(), dname, var+'_'+sample+'.'+ext )
                dest = sample + os.sep + '%s_%s.%s' % (var,i,extNew)
            if verbose:
                print(dest,'-->',src)
            try:
                os.symlink(src,dest)
            except OSError:
                pass

print('sample names: ',sampleNames)
print('field names: ',varNames)
print('time steps: ',len(timesteps),'[',timesteps[indices[0]],'...',timesteps[indices[-1]],']')

