#!/usr/bin/env python
#
# Scrape parameters from specified SOWFA include file
#
import os

def readVarDefinition(line):
    if '//' in line:
        line = line.split('//')
        comment = line[1].strip()
        part1 = line[0].strip()
    else:
        comment = ''
        part1 = line.strip()
    if not part1.endswith(';'):
        return None

    nameval = part1[:-1].split()
    if not len(nameval) == 2:
        return None
    name = nameval[0]
    val = ' '.join(nameval[1:])
    try: 
        val = float(val)
    except ValueError:
        pass

    return (name,val,comment)

def read(fname='setUp',verbose=False):
    params = dict()
    comments = dict()
    if not os.path.isfile(fname):
        print('File not found:',fname)
        return params

    with open(fname,'r') as f:
        processLine = True
        for line in f:
            if line.startswith('//') or line.startswith('#') or line.strip()=='':
                #print('SKIP '+line.rstrip())
                continue
            elif line.startswith('/*'):
                #print('SKIP '+line.rstrip())
                processLine = False
            elif not processLine and '*/' in line:
                #print('SKIP '+line.rstrip())
                processLine = True
            elif processLine:
                p = readVarDefinition(line)
                if not p is None:
                    if verbose: print(p)
                    params[p[0]] = p[1]
                    comments[p[0]] = p[2]

    return params, comments

