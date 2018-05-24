#
# Roughness mapping utilities
#
# Written by Eliot Quon (eliot.quon@nrel.gov), 2018-05-24
#
from __future__ import print_function

def next_name(prefix,newdict):
    names = [ name for name in newdict.keys()
              if name.startswith(prefix) or name.startswith(prefix + ' (')  ]
    return '{:s} ({:d})'.format(prefix,len(names))

def read_landuse_table(tblfile,verbose=True):
    from collections import OrderedDict
    landuse = dict()
    with open(tblfile,'r') as f:
        def read_landuse_def():
            header = f.readline().split(',')
            Ntype = int(header[0])
            Nseason = int(header[1])
            newdict = dict()
            for _ in range(Nseason):
                season = f.readline().strip()
                newdict[season] = OrderedDict()
                for _ in range(Ntype):
                    line = f.readline().strip().split(',')
                    name = line[-1].strip("'")
                    z0 = float(line[4])
                    if name in newdict[season].keys():
                        name = next_name(name, newdict[season])
                    newdict[season][name] = z0
            return newdict
        name = f.readline()
        while not name == '':
            name = name.strip()
            if verbose: print('Read definition:',name)
            landuse[name] = read_landuse_def()
            name = f.readline()
    return landuse

#def read_landuse_table(tblfile,verbose=True):
#    import pandas as pd
#    landuse = []
#    with open(tblfile,'r') as f:
#        def read_landuse_def():
#            line = f.readline().split(',')
#            Ntype = int(line[0])
#            Nseason = int(line[1])
#            header = ['id'] + line[2].strip().strip("'").split() + ['type']
#            dflist = []
#            for _ in range(Nseason):
#                season = f.readline().strip()
#                df = pd.read_csv(f,header=None,names=header,nrows=Ntype)
#                df['season'] = season
#                dflist.append(df)
#                print(df)
#            return pd.concat(dflist)
#        name = f.readline()
#        while not name == '':
#            name = name.strip()
#            if verbose: print('Read definition:',name)
#            df = read_landuse_def()
#            df['definition'] = name
#            landuse.append(df)
#            name = f.readline()
#    return landuse

