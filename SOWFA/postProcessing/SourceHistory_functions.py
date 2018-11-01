#
# Utilities to read in precursor source histories.
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import os
import numpy as np

all_sources = ['UX','UY','UZ','T']

def read_source(source,dpath='./SourceHistory/0'):
    """Read named source located in dpath

    returned S.shape == (Nt,Nz)
    """
    t = []
    dt = []
    S = []
    fpath = os.path.join(dpath,'Source'+source+'History')
    with open(fpath,'r') as f:
        hline = f.readline()
        assert(hline.startswith('Heights (m)'))
        hline = hline.split()
        z = [ float(val) for val in hline[2:] ] # throw out "Heights (m)"
        z = np.array(z)
        f.readline() # throw out header "Time (s) ..."
        for line in f:
            vals = [ float(val) for val in line.split() ]
            t.append( vals[0] )
            dt.append( vals[1] )
            Svals = vals[2:]
            if len(Svals) > 1:
                assert(len(Svals)==len(z))
            S.append(Svals)
        t = np.array(t)
        dt = np.array(dt)
        S = np.array(S)
    return z,t,dt,S


def read_all_sources(dpath):
    S = dict()
    zsave = None
    tsave = None
    dtsave = None
    for Sname in all_sources:
        z,t,dt,S[Sname] = read_source(Sname,dpath=dpath)
        if zsave is None:
            zsave = z
        else:
            assert(np.all(zsave == z))
        if tsave is None:
            tsave = t
        else:
            assert(np.all(tsave == t))
        if dtsave is None:
            dtsave = dt
        else:
            assert(np.all(dtsave == dt))
    return z,t,dt,S

