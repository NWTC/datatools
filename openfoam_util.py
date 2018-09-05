#
# Utilities for handling openfoam files
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
import numpy as np

def read_general(f,line=None,verbose=False):
    """Read scalar, list, or list of lists from currently open file and
    cast resulting lists as arrays where appropriate
    """
    # get a string with a complete definition in one line
    if line is None:
        line = f.readline()
    while (not line.rstrip().endswith(';')) \
            or (not line.count('(') == line.count(')')):
        line += f.readline()
    # clean up the string
    line = line.strip()
    line = line.rstrip(';')
    line = line.replace('(',' ( ') # guarantee parentheses are separated out
    line = line.replace(')',' ) ') # guarantee parentheses are separated out
    line = line.replace('\n',' ')
    line = line.replace('\t',' ')
    # split up the name and data components
    line = list(filter(None, line.split(' '))) # filter removes empty strings
    name = line[0]
    data = []
    if verbose:
        print('key:',name)
        print('raw data (split):',line[1:])
    # at this point, 1-D data looks something like this:
    #   ['(', '2561.5', ')']
    # and 2-D data looks something like this:
    #   ['(', '(', '0.0', '5.26897', ')', '(', '90000.0', '5.26897', ')', ')']
    data = parse_list(line[1:])
    # at this point 1-D data looks something like this:
    #   [2561.5]
    # and 2-D data looks something like this:
    #   [[0.0, 5.26897], [90000.0, 5.26897]]
    if hasattr(data,'__iter__'):
        data = np.array(data)
    return name, data

def parse_list(L,cast=float,verbose=False):
    """Strip out '(' and ')' and replace with list(s)"""
    # single value--do we need to handle this?
    if len(L) == 1:
        L = L[0]
        if cast is not None:
            L = cast(L)
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
        if verbose: print(L)
        # find floats to combine into a list
        for i in range(len(L)):
            if isinstance(L[i], list) or L[i] in ['(',')']:
                continue
            try:
                float(L[i])
            except ValueError:
                continue
            else:
                break
        iend = L.index(')',i)
        # pull out numerical values
        newdata = L[i:iend]
        if cast is not None:
            newdata = [ cast(val) for val in newdata ]
        if verbose: print('new data',newdata)
        # get rid of old list items
        for _ in range(i,iend+1):
            L.pop(i)
        # replace '(' with the new list
        L[i-1] = newdata
        if verbose: print('current list',L)

    return L

def read_all_tables(fname,verbose=True):
    """Read all N-D arrays from the specified file, which is assumed
    to be a basic included input file. Results are returned in a
    dictionary.
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
            elif not line.strip() == '':
                # parse line(s)
                key, val = read_general(f,line)
                data[key] = val
                if verbose: print('read',key)
            # read next line 
            line = f.readline()
    return data

