"""
Class for reading in 'probes' type OpenFOAM sampling

written by Eliot Quon (eliot.quon@nrel.gov)

"""
from __future__ import print_function
import os
import numpy as np

class Probe(object):
    """Stores a time array (t), and field arrays as attributes. The
    fields have shape:
        (N, Nt[, Nd])
    where N is the number of probes and Nt is the number of samples.
    Vectors have an additional dimension to denote vector components.
    Symmetric tensors have an additional dimension to denote tensor components (xx, xy, xz, yy, yz, zz).

    Sample usage:

        from SOWFA.postProcessing.probes import Probe

        # read all probes
        probe = Probe('postProcessing/probe1/0')

        # read specified probes only
        probe = Probe('postProcessing/probe1/0',fields=['U','T'])

        probe.to_csv('probe1.csv')

    """
    def __init__(self,dpath,fields=None):
        """'dpath' is a directory that contains one or more probe files.
        If 'fields' are not explicitly specified, then all probe files
        from the specified directory will be read.
        """
        if fields is None:
            fields = [ fname for fname in os.listdir(dpath)
                        if os.path.isfile(os.path.join(dpath,fname)) ]
        self.fields = []
        for field in fields:
            fpath = os.path.join(dpath,field)
            with open(fpath) as f:
                try:
                    self._read_probe_positions(f)
                except IOError:
                    print('unable to read '+fpath)
                else:
                    self._read_data(f,field)
                    self.fields.append(field)

    def _read_probe_positions(self,f):
        self.pos = []
        line = f.readline()
        while '(' in line and ')' in line:
            line = line.strip()
            assert(line[0]=='#')
            assert(line[-1]==')')
            iprobe = int(line.split()[2])
            i = line.find('(')
            pt = [ float(val) for val in line[i+1:-1].split() ]
            self.pos.append( np.array(pt) )
            line = f.readline()
        if len(self.fields) > 0:
            assert(self.N == len(self.pos))
        else:
            self.N = len(self.pos)
            assert(self.N == iprobe+1)
        self.pos = np.array(self.pos)

    def _read_data(self,f,varname='U'):
        line = f.readline()
        assert(line.split()[1] == 'Time')
        t = []
        F = []
        for line in f:
            line = [ float(val) for val in
                    line.replace('(','').replace(')','').split() ]
            t.append(line[0])
            arr = np.array(line[1:])
            if len(line) == self.N+1:
                # scalar
                F.append(arr)
            elif len(line) == 3*self.N+1:
                # vector
                F.append(arr.reshape((self.N,3),order='C'))
            elif len(line) == 6*self.N+1:
                # symmetric tensor
                F.append(arr.reshape((self.N,6),order='C'))
            else:
                raise IndexError('Unrecognized number of values')
        setattr(self, varname, np.array(F).swapaxes(0,1))
        if len(self.fields) > 0:
            assert(self.Nt == len(self.t))
            assert(np.all(self.t == t))
        else:
            self.t = np.array(t)
            self.Nt = len(self.t)

    def __repr__(self):
        s = 'Times read: {:d} {:s}\n'.format(self.Nt,str(self.t))
        s+= 'Fields read:\n'
        for field in self.fields:
            s+= '  {:s} : {:s}\n'.format(field,
                                         str(getattr(self,field).shape))
        return s

    def as_dataframe(self):
        import pandas as pd
        dflist = []
        for iprobe in range(self.N):
            data = dict(t=self.t)
            for field in self.fields:
                F = getattr(self,field)
                if len(F.shape)==3:
                    # vector
                    for i in range(3):
                        data[field+str(i)] = F[iprobe,:,i]
                else:
                    # scalar
                    data[field] = F[iprobe,:]
            df = pd.DataFrame(data=data)
            df['id'] = iprobe
            dflist.append(df)
        return pd.concat(dflist).set_index(['t','id'])

    def to_csv(self,fname):
        self.as_dataframe().to_csv(fname)

