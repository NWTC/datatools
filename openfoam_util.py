#
# Utilities for handling openfoam files
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
import numpy as np

def _read(f,line=None,debug=False):
    """Read scalar, list, or list of lists from currently open file and
    cast resulting lists as arrays where appropriate
    """
    # get a string with a complete definition in one line
    if line is None:
        line = f.readline()
    while (not ';' in line) \
            or (not line.count('(') == line.count(')')):
        newline = f.readline()
        if newline == '':
            break # EOF
        else:
            line += newline
    # clean up the string
    line = line.strip()
    try:
        line = line[:line.index(';')]
    except ValueError:
        # we reached the end of a list and/or the file, but there was no ';'
        pass
    line = line.replace('(',' ( ') # guarantee parentheses are separated out
    line = line.replace(')',' ) ') # guarantee parentheses are separated out
    line = line.replace('\n',' ')
    line = line.replace('\t',' ')
    # split up the name and data components
    line = list(filter(None, line.split(' '))) # filter removes empty strings
    name = line[0]
    data = []
    if debug:
        print('key:',name)
        print('raw data (split):',line[1:])
    # at this point, 1-D data looks something like this:
    #   ['(', '2561.5', ')']
    # and 2-D data looks something like this:
    #   ['(', '(', '0.0', '5.26897', ')', '(', '90000.0', '5.26897', ')', ')']
    data = of_parse_list(line[1:])
    # at this point 1-D data looks something like this:
    #   [2561.5]
    # and 2-D data looks something like this:
    #   [[0.0, 5.26897], [90000.0, 5.26897]]
    if hasattr(data,'__iter__'):
        data = np.array(data)
    return name, data

def of_parse_list(L,cast=float,debug=False):
    """Strip out '(' and ')' and replace with list(s)"""
    # single value
    if len(L) == 1:
        L = L[0]
        if cast is not None:
            try:
                L = cast(L)
            except ValueError:
                L = L.strip("'").strip('"')
        return L

    # strip first level of list parentheses
    assert((L[0] == '(') and (L[-1] == ')'))
    L = L[1:-1]

    # process 1-D arrays
    if '(' not in L:
        if cast is not None:
            L = [ cast(val) for val in L ]

    # process general arrays
    while '(' in L:
        if debug: print(L)
        # find floats to combine into a list
        for i in range(len(L)):
            if isinstance(L[i], list) or L[i] in ['(',')']:
                continue
            try:
                float(L[i])
            except ValueError:
                # assume arrays/lists only contain numeric items
                continue
            else:
                break
        iend = L.index(')',i)
        # pull out numerical values
        newdata = L[i:iend]
        if cast is not None:
            newdata = [ cast(val) for val in newdata ]
        if debug: print('new data',newdata)
        # get rid of old list items
        for _ in range(i,iend+1):
            L.pop(i)
        # replace '(' with the new list
        L[i-1] = newdata
        if debug: print('current list',L)

    return L

def read_all_defs(fname,verbose=True):
    """Read all definitions, including N-D arrays from the specified
    file, which may be read on the fly during runtime.

    Note: Multi-dimensional OpenFOAM arrays/tables, denoted by nested
    parenthesis, can be read, along with strings and single scalars.
    OpenFOAM dictionaries, denoted by nested curly braces, are NOT
    handled.
    
    Results are returned in a dictionary.
    """
    data = {}
    with open(fname,'r') as f:
        line = f.readline()
        while not line=='': #EOF
            if line.lstrip().startswith('/*'):
                # handle block comments
                while not '*/' in line:
                    line += f.readline()
                if verbose: print(line.rstrip())
            elif line.lstrip().startswith('//'):
                # ignore single line comments
                if verbose: print(line.rstrip())
            elif line.lstrip().startswith('#'):
                # ignore directives
                if verbose: print(line.rstrip())
            elif not line.strip() == '':
                # parse line(s)
                key, val = _read(f,line)
                data[key] = val
                if verbose: print('read',key,val)
            # read next line 
            line = f.readline()
    return data

def of_list(name,L):
    """Return string with a 1-D list/table/array in the OpenFOAM style"""
    s = name + '\n(\n'
    for val in L:
        s += '\t' + str(val) + '\n'
    s += ');\n'
    return s

def of_listlist(name,arr):
    """Return string with a 2-D list/table/array in the OpenFOAM style
    Assume input is a numpy.ndarray
    """
    L = list(arr)
    s = name + '\n(\n'
    for itime,row in enumerate(L):
        s += '\t( '
        s += ' '.join([str(val) for val in row])
        s += ' )\n'
    s += ');\n'
    return s


def seconds_to_datetime(tarray,starttime):
    """Convert an array of times in seconds to an array of datetime
    objects. Start time is a datetime string that can be converted using
    pd.to_datetime().
    """
    import pandas as pd
    tarray = pd.to_timedelta(tarray,unit='s')
    t0 = pd.to_datetime(starttime)
    return tarray + t0


def _get_unique_points_from_list(ylist,zlist,NY=None,NZ=None,order='F'):
    """Detects y and z (1-D arrays) from a list of points on a
    structured grid. Makes no assumptions about the point
    ordering
    """
    ylist = np.array(ylist)
    zlist = np.array(zlist)
    N = len(zlist)
    assert(N == len(ylist))
    if (NY is not None) and (NZ is not None):
        # use specified plane dimensions
        assert(NY*NZ == N)
        y = ylist.reshape((NY,NZ))[:,0]
    elif zlist[1]==zlist[0]:
        # y changes faster, F-ordering
        NY = np.nonzero(zlist > zlist[0])[0][0]
        NZ = int(N / NY)
        assert(NY*NZ == N)
        y = ylist[:NY]
        z = zlist.reshape((NY,NZ),order='F')[0,:]
    elif ylist[1]==ylist[0]:
        # z changes faster, C-ordering
        NZ = np.nonzero(ylist > ylist[0])[0][0]
        NY = int(N / NZ)
        assert(NY*NZ == N)
        z = zlist[:NZ]
        y = ylist.reshape((NY,NZ),order='C')[:,0]
    else:
        print('Unrecognized point distribution')
        print('"y" :',len(ylist),ylist)
        print('"z" :',len(zlist),zlist)
        return ylist,zlist,False
    return y,z,True



