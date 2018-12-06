#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Functions to deal with i/o files. 

VTK structured refers to the following (used in FAST.Farm):

    # <some_header>
    <some_file_descriptor_string>
    ASCII
    DATASET STRUCTURED_POINTS
    DIMENSIONS <nx> <ny> <nz>
    ORIGIN <x0> <y0> <z0>
    SPACING <dx> <dy> <dz>
    POINT_DATA <nx*ny*nz>
    VECTORS <name_of_variable> <type_of_variable_e.g._float>
    
    <value_at_x0> <value_at_y0> <value_at_z0>
    <value_at_x1> <value_at_y0> <value_at_z0>
    <value_at_x2> <value_at_y0> <value_at_z0>    
    
    ...
    
    <value_at_nx-1> <value_at_y0-1> <value_at_z0-1>

"""
from __future__ import print_function
import numpy as np 
import os, glob
import struct

#==============================================================================
# 
#==============================================================================
def _read_vtkStructured_oneFile(vtkPath,verbose=False):
    """
    Reads in a VTK structured file.
    
    Parameters
    ----------
    vtkPath : str,
        full path to VTK file .OR. to a directory full of time sub-directories with vtk files
    verbose : bool,
        whether to print out metadata

    Returns
    -------
    data : list,
        [X,Y,Z,U,V,W] each of which is a 3d array
    meta : dict,
        info about the grid
        
    @author: pdoubraw
    """    
    f = open(vtkPath)
    lines = f.readlines()
    f.close()
    
    nx, ny, nz = [ float(x) for x in lines[4].lstrip().split()[1:] ]
    xo, yo, zo = [ float(x) for x in lines[5].lstrip().split()[1:] ]
    dx, dy, dz = [ float(x) for x in lines[6].lstrip().split()[1:] ]
    npts = int(lines[7].split(" ")[1])
    
    x1d = [xo] if dx==0 else [ xo+i*dx for i in range(int(nx)) ]
    y1d = [yo] if dy==0 else [ yo+i*dy for i in range(int(ny)) ]
    z1d = [ zo+i*dz for i in range(int(nz)) ] # np.arange(zo,zo+dz*nz,dz)           
    
    [Y,X,Z] = np.meshgrid(y1d,x1d,z1d)
    
    U = np.zeros(X.shape)
    V = np.zeros(X.shape)
    W = np.zeros(X.shape)
    
    assert(nx*ny*nz==npts)
    
    # find row index of first numeric value
    for iline,line in enumerate(lines):
        val = lines[iline].lstrip().rstrip().split()[0]    
        try:
            val = float(val)
            if isinstance(val,float):
                row = iline
                break
        except:
            1
            
    # recall that x varies first so this loop must be z->y->x            
    for iz,z in enumerate(z1d):
        for iy,y in enumerate(y1d):
            for ix,x in enumerate(x1d):
                u, v, w = [ float(x) for x in lines[row].lstrip().rstrip().split() ]
                U[ix,iy,iz] = u ; V[ix,iy,iz] = v ; W[ix,iy,iz] = w ; 
                row += 1

    data = [X,Y,Z,U,V,W]                
    meta = {}
    meta['dx'] = dx ; meta['dy'] = dy ; meta['dz'] = dz
    meta['nx'] = nx ; meta['ny'] = ny ; meta['nz'] = nz
    meta['xOrigin'] = xo ; meta['yOrigin'] = yo ; meta['zOrigin'] = zo
    meta['nPts'] = npts

    if verbose:            
        print("dx = {0}".format(dx))
        print("dy = {0}".format(dy))
        print("dz = {0}".format(dz))
        print("nx = {0}".format(nx))
        print("ny = {0}".format(ny))
        print("nz = {0}".format(nz))
        
    return data, meta
#==============================================================================
# 
#==============================================================================
def _read_vtkStructured_manyFiles(vtkPath,t0,dt,nt,verbose=False):
    """
    Read in several *.vtk (VTK Structured) files.
    
    Parameters
    ----------    
    vtkPath : str,
        absolute path to directory where each time subdirectory is (and within each, a vtk structured file)
    t0 : float,
        starting time matching the name of the directory in which the *.xy file will be found
    dt : float,
        time increment in between files/directories
    nt : int,
        number of times to process
    verbose : bool,
        whether to print out messages
        
    Returns
    -------
    data : list,
        [X,Y,Z,U,V,W] each of which is a 3d array
    meta : dict,
        info about the grid
        
    @author: pdoubraw
    """    
    pathNow = os.path.abspath(os.curdir)
    os.chdir(vtkPath)
    times = np.arange(t0,t0+nt*dt,dt)

    for itime,time in enumerate(times):

        timePath        = os.path.abspath(os.path.join(vtkPath,"{0:.3f}".format(time)))
        vtkFile         = glob.glob(os.path.join(timePath,'array*U*.vtk'))[0]
        data, meta      = _read_vtkStructured_oneFile(vtkFile,verbose=False)
        [X,Y,Z,U,V,W]   = data
        
        if verbose:
            print("Reading in {0}...".format(vtkFile))
    
        if itime==0:
            
            x4d = X.copy() ; y4d = Y.copy() ; z4d = Z.copy()
            u4d = U.copy() ; v4d = V.copy() ; w4d = W.copy()
    
            x4d = np.expand_dims(x4d, axis=0)
            y4d = np.expand_dims(y4d, axis=0)
            z4d = np.expand_dims(z4d, axis=0)
            u4d = np.expand_dims(u4d, axis=0)
            v4d = np.expand_dims(v4d, axis=0)
            w4d = np.expand_dims(w4d, axis=0)        
                
        else:
            
            x4d = np.append(x4d,[X],axis=0)
            y4d = np.append(y4d,[Y],axis=0)
            z4d = np.append(z4d,[Z],axis=0)
            u4d = np.append(u4d,[U],axis=0)
            v4d = np.append(v4d,[V],axis=0)
            w4d = np.append(w4d,[W],axis=0)
        
    os.chdir(pathNow)

    return [x4d,y4d,z4d,u4d,v4d,w4d], meta
#==============================================================================
# 
#==============================================================================
def read_vtkStructured(vtkPath,verbose=False):
    """
    Reads in VTK structured data (either one file or a set of files).

    Note: This doesn't work for binary files at the moment --EWQ

    Parameters
    ----------
    vtkPath : str,
        full path to VTK file .OR. to a directory full of time sub-directories with vtk files
    verbose : bool,
        whether to print out metadata

    Returns
    -------
    data : list,
        [X,Y,Z,U,V,W] each of which is a 3d array
    meta : dict,
        info about the grid
        
    @author: pdoubraw    
    """
    extension = os.path.splitext(vtkPath)[-1]
    if extension==".vtk":
        data, meta = _read_vtkStructured_oneFile(vtkPath=vtkPath,verbose=verbose)
    else:
        data, meta = _read_vtkStructured_manyFiles(vtkPath=vtkPath,verbose=verbose)
    return data, meta
#==============================================================================
# 
#==============================================================================
def write_vtkStructured(data,meta,fileOutPath,descStr="PLACEHOLDER",verbose=False):
    """
    Writes data in vtk structured format.
    
    Parameters
    ----------
    data : list,
        [X,Y,Z,U,V,W] each of which is a 3d array
    meta : dict,
        info about the grid
    fileOutPath : str,
        absolute path to vtk file you want to write      
    descStr : str,
        some header string describing what these data are
    """      

    with open(fileOutPath, 'w') as f:      
        f.write('# vtk DataFile Version 3.0\n')  
        f.write('{0}\n'.format(descStr))  
        f.write('ASCII\n')  
        f.write('DATASET STRUCTURED_POINTS\n')  
        f.write('DIMENSIONS {0:d} {1:d} {2:d}\n'.format(int(meta['nx']),
                int(meta['ny']),int(meta['nz'])))  
        f.write('ORIGIN {0:.1f} {1:.1f} {2:.1f}\n'.format(meta['xOrigin'],meta['yOrigin'],meta['zOrigin']))  
        f.write('SPACING {0:.1f} {1:.1f} {2:.1f}\n'.format(meta['dx'],meta['dy'],meta['dz']))  
        f.write('POINT_DATA {0:d}\n'.format(meta['nPts']))  
        f.write('VECTORS vAmb float\n')  
    
        [X,Y,Z,U,V,W]   = data    
        U = np.ravel(U, order='F') ; V = np.ravel(V, order='F') ; W = np.ravel(W, order='F')        
        data = np.zeros((len(U),3))
        data[:,0] = U ; data[:,1] = V ; data[:,2] = W     
        np.savetxt(f,data)

        if verbose:
            print("Saved data to {0}".format(fileOutPath))
            
    return
#==============================================================================
# 
#==============================================================================
def vtk_write_structured_points( f, nx,ny,nz, data,
        datatype=['vector'],
        ds=None,dx=None,dy=None,dz=None,
        origin=(0.0,0.0,0.0),
        dataname=[],
        indexorder='ijk',
        vtk_header='# vtk DataFile Version 2.0',
        vtk_datatype='float',
        vtk_description='really cool data'
    ):
    """Write a VTK dataset with regular topology to file handle 'f'
    written by Eliot Quon (eliot.quon@nrel.gov)
    Note: This should be merged with the existing write_vtkStructured (EWQ)

    Inputs are written with x increasing fastest, then y, then z.

    Example: Writing out two vector fields in one VTK file.
        vtk_write_structured_points(f,nx,ny,nz,[U,V,W,up,vp,wp],ds=1.0,dataname=['mean','fluctuation'])

    Parameters
    ----------
    nx, ny, nz : int
        Data dimensions
    data : list of numpy.ndarray
        The length of this list should correspond to the total number of
        scalars and vector components
    datatype : list
        Acceptable types are 'vector' or 'scalar', and dictate which
        input data correspond to which field
    ds : float, optional
        Default grid spacing; dx,dy,dz may be specified to override
    dx, dy, dz : float, optional
        Specific grid spacings; if ds is not specified, then all three
        must be specified
    origin : list-like, optional
        Origin of the grid
    dataname : list
        List of names for each vector or scalar field in data
    indexorder: str
        Specify the indexing convention (standard: 'ijk', TTUDD: 'jik')

    @author: ewquon
    """
    # calculate grid spacings if needed
    if ds:
        if not dx: dx = ds
        if not dy: dy = ds
        if not dz: dz = ds
    else:
        assert( dx > 0 and dy > 0 and dz > 0 ) 

    # replace shorthand names
    if type(dataname)==str: dataname = [dataname]
    Nvector = 0
    Nscalar = 0
    Nvalues = 0
    for i,name in enumerate(datatype):
        if name[0].lower() == 'v':
            datatype[i] = 'vector'
            Nvector += 1
            Nvalues += 3
        elif name[0].lower() == 's':
            datatype[i] = 'scalar'
            Nscalar += 1
            Nvalues += 1
        else:
            print('unrecognized data type',name)

    # sanity checks
    assert( len(data) == Nvalues )

    # write header
    if 'b' in f.mode:
        binary = True
        import struct
        if bytes is str:
            # python 2
            def b(s):
                return str(s)
        else:
            # python 3
            def b(s):
                return bytes(s,'utf-8')
        f.write(b(vtk_header+'\n'))
        f.write(b(vtk_description+'\n'))
        f.write(b('BINARY\n'))
        f.write(b('DATASET STRUCTURED_POINTS\n'))

        # write out mesh descriptors
        f.write(b('DIMENSIONS {:d} {:d} {:d}\n'.format(nx,ny,nz)))
        f.write(b('ORIGIN {:f} {:f} {:f}\n'.format(origin[0],origin[1],origin[2])))
        f.write(b('SPACING {:f} {:f} {:f}\n'.format(dx,dy,dz)))

        # write out data
        f.write(b('POINT_DATA {:d}\n'.format(nx*ny*nz)))

    else:
        binary = False
        f.write(vtk_header+'\n')
        f.write(vtk_description+'\n')
        f.write('ASCII\n')
        f.write('DATASET STRUCTURED_POINTS\n')

        # write out mesh descriptors
        f.write('DIMENSIONS {:d} {:d} {:d}\n'.format(nx,ny,nz))
        f.write('ORIGIN {:f} {:f} {:f}\n'.format(origin[0],origin[1],origin[2]))
        f.write('SPACING {:f} {:f} {:f}\n'.format(dx,dy,dz))

        # write out data
        f.write('POINT_DATA {:d}\n'.format(nx*ny*nz))

    idx = 0 # data list index
    for idata,outputtype in enumerate(datatype):

        if outputtype=='vector':
            u,v,w = data[idx], data[idx+1], data[idx+2]
            idx += 3
        elif outputtype=='scalar':
            u = data[idx]
            idx += 1
        else: continue

        try:
            #name = dataname[idata]
            name = dataname[idata].replace(' ','_')
        except IndexError:
            name = outputtype+str(idata)

        if outputtype=='vector':
            mapping = { 'i': range(nx), 'j': range(ny), 'k': range(nz) }
            ijkranges = [ mapping[ijk] for ijk in indexorder ]
            if binary:
                f.write(b('{:s}S {:s} {:s}\n'.format(outputtype.upper(),name,vtk_datatype)))
                for k in ijkranges[2]:
                    for j in ijkranges[1]:
                        for i in ijkranges[0]:
                            f.write(struct.pack('>fff', u[i,j,k], v[i,j,k], w[i,j,k])) # big endian
            else: #ascii
                f.write('{:s}S {:s} {:s}\n'.format(outputtype.upper(),name,vtk_datatype))
                for k in ijkranges[2]:
                    for j in ijkranges[1]:
                        for i in ijkranges[0]:
                            f.write(' {:f} {:f} {:f}\n'.format(u[i,j,k], v[i,j,k], w[i,j,k]))
        elif outputtype=='scalar':
            if binary:
                f.write(b('{:s}S {:s} {:s}\n'.format(outputtype.upper(),name,vtk_datatype)))
                f.write(b('LOOKUP_TABLE default\n'))
                for k in range(nz):
                    for j in range(ny):
                        for i in range(nx):
                            #f.write(struct.pack('f',u[j,i,k])) # native endianness
                            f.write(struct.pack('>f',u[j,i,k])) # big endian
            else:
                f.write('{:s}S {:s} {:s}\n'.format(outputtype.upper(),name,vtk_datatype))
                f.write('LOOKUP_TABLE default\n')
                for k in range(nz):
                    for j in range(ny):
                        for i in range(nx):
                            f.write(' {:f}\n'.format(u[j,i,k]))


def vtk_read_binary_structured_points(fname,dtype=np.float32,verbose=True):
    """Read VTK dataset written with vtk_write_structured_points
    Note: At the moment read_vtkStructured doesn't properly handle
    binary files
    """
    if verbose:
        def readecho(): print(f.readline().strip())
    else:
        def readecho(): f.readline()
    prec = np.dtype(dtype).itemsize
    vectorData = dict()
    scalarData = dict()
    with open(fname,'rb') as f:
        readecho() # header
        readecho() # description
        readecho() # file mode
        readecho() # expected: DATASET
        dims = [int(val) for val in f.readline().split()[1:]]
        origin = [float(val) for val in f.readline().split()[1:]]
        spacing = [float(val) for val in f.readline().split()[1:]]
        N = np.prod(dims)
        readecho() # expected: POINT_DATA
        newdataline = f.readline()
        while not newdataline=='':
            fieldtype, name, datatype = newdataline.split()
            print('Processing {}field {} (dtype={})'.format(fieldtype.lower().strip('s'),name,datatype))
            if fieldtype.lower().startswith('vector'):
                #vectorData[name] = np.zeros([3]+dims,dtype=dtype)
                data = struct.unpack('>{:d}f'.format(3*N),f.read(3*N*prec))
                vectorData[name] = np.array(data,dtype=dtype).reshape([3]+dims,order='F')
            elif fieldtype.lower().startswith('scalar'):
                #scalarData[name] = np.zeros(dims,dtype=dtype)    
                scalarData[name] = np.array(data,dtype=dtype).reshape(dims,order='F')
            newdataline = f.readline()
    if verbose:
        print('Read scalar data:',scalarData.keys())
        print('Read vector data:',vectorData.keys())
    meta = dict()
    meta['dx'] = spacing[0]
    meta['dy'] = spacing[1]
    meta['dz'] = spacing[2]
    meta['nx'] = dims[0]
    meta['ny'] = dims[1]
    meta['nz'] = dims[2]
    meta['xOrigin'] = origin[0]
    meta['yOrigin'] = origin[1]
    meta['zOrigin'] = origin[2]
    meta['nPts'] = N
    return scalarData,vectorData,meta
