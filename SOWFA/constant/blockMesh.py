#
# Module for assisting with blockMesh grid generation
#
# Written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import numpy as np

def grow_mesh(N, d0, r):
    """Calculate the cell spacings given a number of cells (N), an 
    initial spacing (d0), and a normal growth rate (r).

    Spacings are calculated as follows:
        d_i = d_0 * r^i,  for i in 0 .. N-1
    """
    spacings = d0 * np.array([r**i for i in range(N)])
    print('spacings : {:g} .. {:g}'.format(d0,spacings[-1]))
    return spacings

def points_from(spacings, x0=0.0):
    x = x0 * np.ones((len(spacings)+1,))
    x[1:] += np.cumsum(spacings)
    return x

def spacings(L, N, ratio=1):
    """Calculate the blockmesh spacings that would result given a length
    (L), number of cells (N), and simpleGrading ratio.
    """
    r = np.exp(np.log(ratio) / (N-1))
    print('growth rate = ',r)
    d0 = L / np.sum([r**i for i in range(N)])
    return grow_mesh(N, d0, r)

def start_end(L, N, ratio):
    r = np.exp(np.log(ratio) / (N-1))
    d0 = L / np.sum([r**i for i in range(N)])
    return d0, ratio*d0

def estimate(d, d0, r, L=None, round_to=10., output='blockMesh'):
    """Estimate the blockmesh inputs that would give the desired final
    point spacing (d) over a specified distance (L) with approximately
    the normal growth rate (r). The actual growth rate will not exactly
    match the input unless d_n == d0 * r^(N-1). 
    
    The initial spacing (d0) is used to estimate the number of points
    within the region, and is then adjusted to give the desired final
    spacing (d) over distance (L).

    If L is not specified, then a reasonable estimate of blockMesh
    inputs is provided with L rounded to a "nice" number.

    Returns estimated N, d0, r
    """
    # first estimate the number of cells to get closed to the initial
    # and final spacings given a particular growth rate
    N_est = np.log(d/d0) / np.log(r) + 1
    N = int(np.ceil(N_est))
    print('estimated number of cells is {:g} ~= {:d}'.format(N_est, N))
    L1 = np.sum(grow_mesh(N,d0,r))
    print('resultant layer height : ',L1)
    # we assume that the estimated number of cells is pretty close to
    # what it should be and leave this fixed; since the growth rate
    # r = f(d0,d1,N), we should update it assuming the other parameters
    # are also fixed
    r_approx = np.exp(np.log(d/d0) / (N-1))
    print('actual growth rate: ',r_approx)
    # we're basically done at this point, but the resulting distance 
    # spanned by the points probably isn't a nice floating point number
    if L is None:
        L = np.round(L1/round_to) * round_to
    print('calculating d0 for L =',L)
    d0_approx = L / np.sum([r_approx**i for i in range(N)])
    print('adjusted initial spacing :',d0_approx)
    if output=='growth':
        return N, d0_approx, r_approx
    else:
        return L, N, r_approx**(N-1)

