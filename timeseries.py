#!usr/bin/env python
from __future__ import print_function
import os

def pretty_list(strlist,indent=2,sep='\t',width=80):
    """For formatting long lists of strings of arbitrary length
    """
    sep = sep.expandtabs()
    max_item_len = max([len(s) for s in strlist])
    items_per_line = int((width - (indent+max_item_len)) / (len(sep)+max_item_len) + 1)
    Nlines = int(len(strlist) / items_per_line)
    extraline = (len(strlist) % items_per_line) > 0
    fmtstr = '{{:{:d}s}}'.format(max_item_len)
    strlist = [ fmtstr.format(s) for s in strlist ] # pad strings so that they're all the same length
    finalline = ''
    for line in range(Nlines):
        ist = line*items_per_line
        finalline += indent*' ' + sep.join(strlist[ist:ist+items_per_line]) + '\n'
    if extraline:
        finalline += indent*' ' + sep.join(strlist[Nlines*items_per_line:]) + '\n'
    return finalline


class TimeSeries(object):
    """Object for holding general time series data in a single
    directory

    Written by Eliot Quon (eliot.quon@nrel.gov)

    Sample usage:
        from datatools.timeseries import ATimeSeries
        ts = TimeSeries('/path/to/data',prefix='foo',suffix='.bar')
        for fname in ts:
            do_something(fname)
        for t,fname in ts.itertimes():
            do_something(t,fname)
    """

    def __init__(self,
                 datadir='.',
                 prefix=None, suffix='',
                 dt=1.0, t0=0.0,
                 verbose=False):
        """Collect data from specified directory, for files with a
        given prefix and optional suffix. For series with integer time
        step numbers, dt can be specified (with t0 offset) to determine
        the time for each snapshot.
        """
        # default initialization for inherited objects
        self.datadir = os.path.abspath(datadir)
        self.filelist = None
        self.times = []
        self.verbose = verbose

        if prefix is not None:
            # handle standard time series in a single directory
            self.filelist = []
            self.times = []
            self.dt = dt
            self.t0 = t0
            for f in os.listdir(self.datadir):
                if (os.path.isfile(os.path.join(self.datadir,f))) \
                        and f.startswith(prefix) \
                        and f.endswith(suffix):
                    fpath = os.path.join(self.datadir,f)
                    self.filelist.append(fpath)
                    val = f[len(prefix):]
                    if len(suffix) > 0:
                        val = val[:-len(suffix)]
                    try:
                        self.times.append(t0 + dt*float(val))
                    except ValueError:
                        print('Prefix and/or suffix are improperly specified')
                        print('  attempting to cast value: '+val)
                        print('  for file: '+fpath)
                        break
            self.Ntimes = len(self.filelist)
            if self.Ntimes == 0:
                print('Warning: no matching files were found')

            # sort by output time
            iorder = [kv[0] for kv in sorted(enumerate(self.times),key=lambda x:x[1])]
            self.filelist = [self.filelist[i] for i in iorder]
            self.times = [self.times[i] for i in iorder]

    def __len__(self):
        return len(self.filelist)

    def __getitem__(self,i):
        return self.filelist[i]

    def __iter__(self):
        self.lastfile = -1  # reset iterator index
        return self

    def __next__(self):
        if self.filelist is None:
            raise StopIteration('file list is empty')
        self.lastfile += 1
        if self.lastfile >= self.Ntimes:
            raise StopIteration
        else:
            return self.filelist[self.lastfile]

    def next(self):
        # for Python 2 compatibility
        return self.__next__()
            
    def itertimes(self):
        return zip(self.times,self.filelist)


class SOWFATimeSeries(TimeSeries):
    """Object for holding general time series data stored in multiple
    time subdirectories, e.g., as done in OpenFOAM.

    Written by Eliot Quon (eliot.quon@nrel.gov)

    Sample usage:
        from datatools.timeseries import SOWFATimeSeries
        ts = SOWFATimeSeries('/path/to/data',filename='U')
    """

    def __init__(self,datadir='.',filename=None,verbose=True):
        """Collect data from subdirectories, assuming that subdirs
        have a name that can be cast as a float
        """
        super(self.__class__,self).__init__(datadir=datadir,verbose=verbose)
        self.dirlist = []
        self.timenames = []
        self.filename = filename

        # process all subdirectories
        subdirs = [ os.path.join(self.datadir,d)
                    for d in os.listdir(self.datadir)
                    if os.path.isdir(os.path.join(self.datadir,d)) ]
        for path in subdirs:
            dname = os.path.split(path)[-1]
            try:
                tval = float(dname)
            except ValueError:
                continue
            self.times.append(tval)
            self.dirlist.append(path)
        self.Ntimes = len(self.dirlist)
    
        # sort by output time
        iorder = [kv[0] for kv in sorted(enumerate(self.times),key=lambda x:x[1])]
        self.dirlist = [self.dirlist[i] for i in iorder]
        self.times = [self.times[i] for i in iorder]

        # check that all subdirectories contain the same files
        self.timenames = os.listdir(self.dirlist[0])
        for d in self.dirlist:
            if not os.listdir(d) == self.timenames:
                print('Warning: not all subdirectories contain the same files')
                break
        if verbose:
            self.outputs() # print available outputs

        # set up file list
        if filename is not None:
            self.update_filelist(filename)

    def update_filelist(self,filename):
        """Update file list for iteration"""
        self.filelist = []
        for path in self.dirlist:
            fpath = os.path.join(path,filename)
            if os.path.isfile(fpath):
                self.filelist.append(fpath)
            else:
                raise IOError(fpath+' not found')

    def outputs(self,prefix=''):
        """Print available outputs for the given data directory"""
        selected_output_names = [ name for name in self.timenames if name.startswith(prefix) ]
        if self.verbose:
            if prefix:
                print('Files starting with "{}" in each subdirectory:'.format(prefix))
            else:
                print('Files in each subdirectory:')
            #print('\t'.join([ '    '+name for name in selected_output_names ]))
            print(pretty_list(sorted(selected_output_names)))
        return selected_output_names

    def __repr__(self):
        return str(self.Ntimes) + ' time subdirectories located in ' + self.datadir

