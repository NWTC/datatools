#!/usr/bin/env python
import os

def pretty_list(strlist,indent=2,sep='\t',width=80):
    """For formatting long lists of strings of arbitrary length
    """
    sep = sep.expandtabs()
    max_item_len = max([len(s) for s in strlist])
    items_per_line = (width - (indent+max_item_len)) / (len(sep)+max_item_len) + 1
    Nlines = len(strlist) / items_per_line
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
    """ Object for holding general time series data which may be stored
    in multiple time subdirectories

    Written by Eliot Quon (eliot.quon@nrel.gov)

    Sample usage:
        from datatools.timeseries import TimeSeries
        ts = TimeSeries('/path/to/data',filename='U')
        for fname in ts:
            do_something(fname)
    """

    def __init__(self,datadir='.',filename=None,verbose=True):
        """ Collect data from subdirectories, assuming that subdirs
        have a name that can be cast as a float
        """
        self.dataDir = os.path.abspath(datadir)
        self.filename = filename
        self.outputTimes = []
        self.outputNames = []
        self.dirList = []
        self.fileList = None
        self.lastFile = -1  # for iterator
        self.verbose = verbose

        # process all subdirectories
        subdirs = [ os.path.join(self.dataDir,d)
                    for d in os.listdir(self.dataDir)
                    if os.path.isdir(os.path.join(self.dataDir,d)) ]
        for path in subdirs:
            dname = os.path.split(path)[-1]
            try:
                tval = float(dname)
            except ValueError:
                continue
            self.outputTimes.append(tval)
            self.dirList.append(path)
        self.Ntimes = len(self.dirList)
    
        # sort by output time
        iorder = [kv[0] for kv in sorted(enumerate(self.outputTimes),key=lambda x:x[1])]
        self.dirList = [self.dirList[i] for i in iorder]
        self.outputTimes = [self.outputTimes[i] for i in iorder]

        # check that all subdirectories contain the same files
        self.outputNames = os.listdir(self.dirList[0])
        for d in self.dirList:
            if not os.listdir(d) == self.outputNames:
                print 'Warning: not all subdirectories contain the same files'
                break
        if verbose:
            self.outputs() # print available outputs

        # set up file list
        if filename is not None:
            self.setFilename(filename)

    def setFilename(self,filename):
        """Update file list for iteration"""
        self.lastFile = -1  # reset iterator index
        self.fileList = []
        for path in self.dirList:
            fpath = os.path.join(path,filename)
            if os.path.isfile(fpath):
                self.fileList.append(fpath)
            else:
                raise IOError(fpath+' not found')

    def outputs(self,prefix=''):
        """Print available outputs for the given data directory"""
        selectedOutputNames = [ name for name in self.outputNames if name.startswith(prefix) ]
        if self.verbose:
            if prefix:
                print 'Files starting with "'+prefix+'" in each subdirectory:'
            else:
                print 'Files in each subdirectory:'
            #print '\t'.join([ '    '+name for name in selectedOutputNames ])
            print pretty_list(sorted(selectedOutputNames))
        return selectedOutputNames

    def __repr__(self):
        return str(self.Ntimes) + ' time subdirectories located in ' + self.dataDir

    def __len__(self):
        return len(self.fileList)

    def __getitem__(self,i):
        return self.fileList[i]

    def __iter__(self):
        return self

    def next(self):
        if self.fileList is None:
            raise StopIteration('Need to set filename before iterating')
        self.lastFile += 1
        if self.lastFile >= self.Ntimes:
            raise StopIteration
        else:
            return self.fileList[self.lastFile]
            

