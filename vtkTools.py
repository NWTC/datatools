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

import numpy as np 
import os, glob

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
    z1d = np.arange(zo,zo+dz*nz,dz)           
    
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
        print "dx = {0}".format(dx)
        print "dy = {0}".format(dy)
        print "dz = {0}".format(dz)
        print "nx = {0}".format(nx)
        print "ny = {0}".format(ny)
        print "nz = {0}".format(nz)
        
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
            print "Reading in {0}...".format(vtkFile)
    
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
    f = open("top.txt", 'w')
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
    f.close()    
    
    [X,Y,Z,U,V,W]   = data    
    U = np.ravel(U, order='F') ; V = np.ravel(V, order='F') ; W = np.ravel(W, order='F')        
    data = np.zeros((len(U),3))
    data[:,0] = U ; data[:,1] = V ; data[:,2] = W     
    np.savetxt("bot.txt",data)
    os.system("cat top.txt bot.txt > {0}".format(fileOutPath))
    os.remove("top.txt")
    os.remove("bot.txt")    

    if verbose:
        print "Saved data to {0}".format(fileOutPath)
    return