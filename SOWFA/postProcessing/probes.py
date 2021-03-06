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
        probe = Probe('postProcessing/probe1/')

        # read specified probes only
        probe = Probe('postProcessing/probe1/',fields=['U','T'])

        probe.to_csv('probe1.csv')

    """
    def __init__(self,dpath=None,**kwargs):
        """'Find and process all time directories of a particular probe with path dpath"""
        self._processed = []
        self.simTimeDirs = [] #output time names
        self.simStartTimes = [] # start or restart simulation times
        self.imax = None # for truncating time series

        if not dpath:
            dpath = '.'

        if dpath[-1] == os.sep: dpath = dpath[:-1] # strip trailing slash
        
        # find results
        listing = os.listdir(dpath)
        for dirname in listing:
            if not os.path.isdir(dpath+os.sep+dirname): continue
            try:
                startTime = float(dirname)
                self.simTimeDirs.append( dpath+os.sep+dirname )
                self.simStartTimes.append( startTime )
            except ValueError:
                # dirname is not a number
                pass

        # sort results
        self.simTimeDirs = [ x[1] for x in sorted(zip(self.simStartTimes,self.simTimeDirs)) ]
        self.simStartTimes.sort()

        print('Simulation (re)start times:',self.simStartTimes)

        # process all output dirs
        if len(self.simTimeDirs) > 0:
            self._processdirs( self.simTimeDirs, **kwargs )
        else:
            print('No probe time directories found!')
            
        self._trim_series_if_needed(self._processed)


    def _processdirs(self,
                     tdirList,
                     varList=['U','T'],
                     trimOverlap=True
                    ):
        """Reads all files within a probe output time directory.
        An object attribute corresponding to the probe output name
        is updated, e.g.:
            ${timeDir}/U is appended to the array self.U
        """
        if isinstance( varList, (str,) ):
            if varList.lower()=='all':
                # special case: read all vars
                outputs = [ fname for fname in os.listdir(tdirList[0])
                                if os.path.isfile(tdirList[0]+os.sep+fname) ]
            else: # specified single var
                outputs = [varList]
        else: #specified list
            outputs = varList
        
        # process all data
        selected = []
        for field in outputs:
            arrays = [ self._read_probe_file( tdir + os.sep + field ) for tdir in tdirList ]

            # combine into a single array and trim end of time series
            # (because simulations that are still running can have different
            # array lengths)
            newdata = np.concatenate(arrays)[:self.imax,:]

            # get rid of overlapped data for restarts
            if trimOverlap:
                if len(selected) == 0:
                    # create array mask
                    tpart = [ array[:,0] for array in arrays ]
                    for ipart,tcutoff in enumerate(self.simStartTimes[1:]):
                        selectedpart = np.ones(len(tpart[ipart]),dtype=bool)
                        try:
                            iend = np.nonzero(tpart[ipart] >= tcutoff)[0][0]
                        except IndexError:
                            # clean restart
                            pass
                        else:
                            # previous simulation didn't finish; overlapped data
                            selectedpart[iend:] = False 
                        selected.append(selectedpart)
                    # last / currently running part
                    selected.append(np.ones(len(tpart[-1]),dtype=bool))
                    selected = np.concatenate(selected)[:self.imax]
                    assert(len(selected) == len(newdata[:,0]))
                elif not (len(newdata[:,0]) == len(selected)):
                    # if simulation is still running, subsequent newdata may
                    # be longer
                    self.imax = min(len(selected), len(newdata[:,0]))
                    selected = selected[:self.imax]
                    newdata = newdata[:self.imax,:]
                # select only unique data
                newdata = newdata[selected,:]

            # reshape field into (Nt,Nz[,Nd]) and set as attribute
            # - note: first column of 'newdata' is time
            # - note: old behavior was to return (Nz,Nt,Nd) for Nd >= 1
            if newdata.shape[1] == self.N+1:
                # scalar
                setattr( self, field, newdata[:,1:] )
            elif newdata.shape[1] == 3*self.N+1:
                # vector
                setattr( self, field, newdata[:,1:].reshape((newdata.shape[0],self.N,3),order='C') )
            elif newdata.shape[1] == 6*self.N+1:
                # symmetric tensor
                setattr( self, field, newdata[:,1:].reshape((newdata.shape[0],self.N,6),order='C') )
            else:
                raise IndexError('Unrecognized number of values')
            self._processed.append(field)
            print('  read',field)        # set time arrays
            
        self.t = newdata[:,0]
        self.Nt = len(self.t)


    def _read_probe_file(self,fpath):
        with open(fpath) as f:
            try:
                self._read_probe_positions(f)
            except IOError:
                print('unable to read '+fpath)
            else:
                array = self._read_data(f)
        return array


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
        if len(self._processed) > 0: # assert that all fields have same number of probes
            assert(self.N == len(self.pos))
        else: # first field: set number of probes in self.N
            self.N = len(self.pos)
            assert(self.N == iprobe+1)
        self.pos = np.array(self.pos)


    def _read_data(self,f):
        line = f.readline()
        assert(line.split()[1] == 'Time')
        out = []
        for line in f:
            line = [ float(val) for val in
                    line.replace('(','').replace(')','').split() ]
            out.append(line)
        return np.array(out)


    def _trim_series_if_needed(self,fields_to_check=None):
        """check for inconsistent array lengths and trim if needed"""
        if fields_to_check is None:
            fields_to_check = self._processed
        for field in fields_to_check:
            try:
                getattr(self,field)
            except AttributeError:
                print('Skipping time series length check for unknown field: ',
                      field)
                fields_to_check.remove(field)
        field_lengths = [ getattr(self,field).shape[0] for field in fields_to_check ]
        if np.min(field_lengths) < np.max(field_lengths):
            self.imax = np.min(field_lengths)
            # need to prune arrays
            print('Inconsistent averaging field lengths... is simulation still running?')
            print('  truncated field histories from',np.max(field_lengths),'to',self.imax)
            self.t = self.t[:self.imax]
            self.Nt = len(self.t)
            for field in fields_to_check:
                Ndim = len(getattr(self,field).shape)
                if Ndim == 2:
                    # scalar
                    setattr(self, field, getattr(self,field)[:self.imax,:])
                elif Ndim == 3:
                    # vector/tensor
                    setattr(self, field, getattr(self,field)[:self.imax,:,:])
                else:
                    print('Unknown field type ',field)


    def __repr__(self):
        s = 'Times read: {:d} {:s}\n'.format(self.Nt,str(self.t))
        s+= 'Fields read:\n'
        for field in self._processed:
            s+= '  {:s} : {:s}\n'.format(field,
                                         str(getattr(self,field).shape))
        return s


    #============================================================================
    #
    # DATA I/O
    #
    #============================================================================

    def to_pandas(self):
        import pandas as pd
        dflist = []
        for iprobe in range(self.N):
            data = dict(t=self.t)
            for field in self._processed:
                F = getattr(self,field)
                if len(F.shape)==2:
                    # scalar
                    data[field] = F[:,iprobe]
                elif F.shape[2]==3:
                    # vector
                    for i,name in enumerate(['x','y','z']):
                        data[field+name] = F[:,iprobe,i]
                elif F.shape[2]==6:
                    # symmetric tensor
                    for i,name in enumerate(['xx','xy','xz','yy','yz','zz']):
                        data[field+name] = F[:,iprobe,i]
            df = pd.DataFrame(data=data)
            #df['id'] = iprobe
            df['z'] = self.pos[iprobe,2]
            dflist.append(df)
        #return pd.concat(dflist).set_index(['t','id'])
        return pd.concat(dflist).sort_values(['t','z']).set_index(['t','z'])

    def to_csv(self,fname):
        self.to_pandas().to_csv(fname)

    def to_netcdf(self,fname):
        long_names = {'T': 'Potential temperature',
                      'Ux': 'U velocity component',
                      'Uy': 'V velocity component',
                      'Uz': 'W velocity component',
                      }
        units = {'T': 'K',
                 'Ux': 'm s-1',
                 'Uy': 'm s-1',
                 'Uz': 'm s-1',
                }

        print('Dumping data to',fname)
        import netCDF4
        f = netCDF4.Dataset(fname,'w')
        f.createDimension('time',len(self.t))
        f.createDimension('z',self.pos.shape[1])

        times = f.createVariable('time', 'float', ('time',))
        times.long_name = 'Time'
        times.units = 's'
        times[:] = self.t

        heights = f.createVariable('z', 'float', ('z',))
        heights.long_name = 'Height above ground level'
        heights.units = 'm'
        heights[:] = self.pos[:,2]

        for var in self._processed:
            F = getattr(self,var)
            if len(F.shape)==2:
                # scalar
                varnames = [var,]
                F = F[:,:,np.newaxis]
            elif F.shape[2]==3:
                # vector
                varnames = [var+name for name in ['x','y','z']]
            elif F.shape[2]==6:
                # symmetric tensor
                varnames = [var+name for name in ['xx','xy','xz','yy','yz','zz']]
            
            for i, varname in enumerate(varnames):
                field = f.createVariable(varname, 'float', ('time','z'))
                try:
                    field.long_name = long_names[varname]
                except KeyError:
                    # Use var name as description
                    field.long_name = varname
                try:
                    field.units = units[varname]
                except KeyError:
                    # Units unknown
                    pass
                field[:] = F[:,:,i]
        f.close()
