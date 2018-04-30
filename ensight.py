#!/usr/bin/env python
import numpy as np

def read_mesh(fname):
    """Read a single part from an ascii ensight mesh file"""
    with open(fname,'r') as f:
        for _ in range(8): f.readline()
        N = int(f.readline())
        x = np.zeros((N,))
        y = np.zeros((N,))
        z = np.zeros((N,))
        for i in range(N):
            x[i] = float(f.readline())
        for i in range(N):
            y[i] = float(f.readline())
        for i in range(N):
            z[i] = float(f.readline())
    return x,y,z

def read_scalar(fname,N):
    """Read a single part from an ascii ensight soln file"""
    s = np.zeros((N,))
    with open(fname,'r') as f:
        fieldtype = f.readline().strip()
        assert(fieldtype == 'scalar')
        for _ in range(3): f.readline()
        for i in range(N):
            s[i] = float(f.readline())
    return s                       

def read_vector(fname,N):
    """Read a single part from an ascii ensight soln file"""
    u = np.zeros((N,))
    v = np.zeros((N,))
    w = np.zeros((N,))
    with open(fname,'r') as f:
        fieldtype = f.readline().strip()
        assert(fieldtype == 'vector')
        for _ in range(3): f.readline()
        for i in range(N):
            u[i] = float(f.readline())
        for i in range(N):
            v[i] = float(f.readline())
        for i in range(N):
            w[i] = float(f.readline())
    return u,v,w                       

