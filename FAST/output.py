#!/usr/bin/env python
#
# Helper module for postprocessing FAST outputs
# Tested with OpenFAST v8.17.00a-bjj
#
# Written by Eliot Quon (eliot.quon@nrel.gov) 2017-07-31
#
import numpy as np
import matplotlib.pyplot as plt

def read(outputFile,**kwargs):
    return FASToutput(outputFile,**kwargs)

class FASToutput(object):
    # TODO: override dict instead of generic object

    def __init__(self,fname=None,verbose=True):
        # inputs
        self.fname = fname
        self.verbose = verbose
        # initialize members
        self.outputs = None
        self.units = None
        self.output_units = dict()
        self.Noutputs = 0
        self.N = 0
        # read output file
        if fname is not None:
            self._readFASToutput(fname)
            if verbose:
                self.printStats()

    def __getitem__(self, key):
        if key in self.outputs:
            return getattr(self, key)
        else:
            raise KeyError('Requested key \'{:s}\' not in {}'.format(key,self.outputs))

    def _readFASToutput(self,fname,Nheaderlines=6):
        if self.verbose: print 'Reading header info from',fname
        with open(fname,'r') as f:
            if self.verbose:
                for _ in range(Nheaderlines): print f.readline().strip()
            else:
                for _ in range(Nheaderlines): f.readline()
            # read names of output quantities
            self.outputs = f.readline().split()
            # read units for each output quantity
            self.units = [ s.strip('()') for s in f.readline().split() ]
            self.Noutputs = len(self.outputs)
            assert(self.Noutputs == len(self.units))
            for iline,_ in enumerate(f): pass
        self.N = iline + 1
        # read data
        if self.verbose: print 'Reading data...'
        data = np.loadtxt(fname,skiprows=Nheaderlines+2)
        for i,output in enumerate(self.outputs):
            setattr(self,output,data[:,i])
            self.output_units[output] = self.units[i]
        assert(len(self.Time) == self.N)
        # aliases
        self._setAlias('t','Time')
        self._setAlias('P','RotPwr')
        self._setAlias('T','RotThrust','LSShftFxa','LSShftFxs','LSSGagFxa','LSSGagFxs')
        self._setAlias('rpm','LSSTipVxa')
        self._setAlias('genspd','HSShftV')
        self._setAlias('pitch1','PtchPMzc1')
        self._setAlias('pitch2','PtchPMzc2')
        self._setAlias('pitch3','PtchPMzc3')
        self._setAlias('pitch','pitch1')

    def addOutput(self,name,data,units=None):
        if name not in self.outputs:
            self.outputs.append(name)
            setattr(self,name,data)
            if units is not None: self.output_units[name] = units
        else:
            print 'Output',name,'already exists'

    def _setAlias(self,name,*aliases):
        for alias in aliases:
            try:
                data = getattr(self,alias)
            except AttributeError:
                continue
            else:
                self.addOutput(name, data)
                if self.verbose: print '  set',name,'-->',alias
                return
        if self.verbose: print 'Outputs for alias',name,'do not exist:',aliases

    def printStats(self):
        print 'Output       Units    Mean         Min          Max          Stdev'
        print '------------ -------- ------------ ------------ ------------ ------------'
        for output,units in zip(self.outputs,self.units):
            data = getattr(self,output)
            print '{:12s} {:8s} {:12g} {:12g} {:12g} {:12g}'.format(
                        output,
                        units,
                        np.min(data),
                        np.max(data),
                        np.mean(data),
                        np.std(data)
                    )
        print ''

    def plot(self,outputName,*args,**kwargs):
        data = getattr(self,outputName)
        Ndata = len(data)
        Nt = len(self.t)
        if Ndata == Nt:
            time = self.t
        elif Ndata < Nt:
            Navg = Nt - Ndata
            time = self.t[Navg/2:-Navg/2]
        plt.plot(time, data, *args, **kwargs)


    def running_mean(self,outputName,Tavg):
        N = int(Tavg / (self.t[1]-self.t[0]))
        data = getattr(self,outputName)
        mean = np.convolve(data, np.ones((N,))/N, mode='valid')
        newname = outputName + '_mean'
        self.addOutput(newname,mean)
        print 'Averaged {:s} with {:f} s window (N={:d})'.format(outputName,self.t[N+1]-self.t[0],N)
        return mean


    def low_pass_filtered_mean(self,outputName,fc=np.inf,order=2):
        # Example: https://gist.github.com/junzis/e06eca03747fc194e322
        from scipy.signal import butter, lfilter
        data = getattr(self,outputName)
        fs = 1./(self.t[1] - self.t[0])
        cutoff_norm = fc / (0.5*fs) # Wn normalized from 0 to 1, where 1 is the Nyquist frequency
        b,a = butter(order, cutoff_norm, btype='lowpass', analog=False, output='ba')
        filtered_data = lfilter(b, a, data)
        newname = outputName + '_mean'
        self.addOutput(newname,filtered_data)
        print 'Filtered {:s} with cutoff freq={:f} Hz, order={:d}'.format(outputName,fc,order)
        return filtered_data


    def fluctuations(self,outputName,meanName=None):
        if meanName is None:
            meanName = outputName + '_mean'
            if not hasattr(self,meanName):
                self.low_pass_filtered_mean(outputName)
        data = getattr(self,outputName)
        mean = getattr(self,meanName)
        if len(data) > len(mean):
            Navg = len(data) - len(mean)
            data = data[Navg/2:-Navg/2]
        fluc = data - mean
        newname = outputName + '_fluc'
        self.addOutput(newname,fluc)
        return fluc


    def vector_magnitude(self,outputname,*components):
        mag_sq = 0.0
        for comp in components:
            mag_sq += self[comp]**2
        units = self.output_units[components[0]]
        self.addOutput(outputname, mag_sq**0.5, units=units) 

