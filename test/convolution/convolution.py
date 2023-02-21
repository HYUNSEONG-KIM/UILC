from typing import Tuple, Literal, Callable

import numpy as np
from numpy.lib.stride_tricks import as_strided

def toeplitz(c, r):
    # Imgrated source code from scipy `toeplitz` function
    # BSD-licensed
    c = np.asarray(c).ravel()
    if r is None:
        r = c.conjugate()
    else:
        r = np.asarray(r).ravel()
    # Form a 1-D array containing a reversed c followed by r[1:] that could be
    # strided to give us toeplitz matrix.
    vals = np.concatenate((c[::-1], r[1:]))
    out_shp = len(c), len(r)
    n = vals.strides[0]
    return as_strided(vals[len(c)-1:], shape=out_shp, strides=(-n, n)).copy()



def _expand_matrix(data, dim_kernel, dim_crop, edge_param):
    n, m = data.shape
    l, k = dim_kernel
    cr, cc = dim_crop

    e_method_name, e_method_param = edge_param

    er = l - cr
    ec = k -cc

    return None

def convolution2d(  
    data:np.ndarray, 
    filter:np.ndarray,
    boundary_cropping:Tuple[int, int]=(1, 1),
    edge_handling:Literal["extend", "wrap", "mirror","constant"]="constant",
    constant = 0
    ):

    data.shape

    return None

def convert_to_toeplitz_system(data, filter, 
    boundary_cropping:Tuple[int, int]=(1, 1),
    edge_handling:Literal["extend", "wrap", "mirror","constant"]="constant",
    constant = 0, preserve_filter=False) -> Tuple[np.ndarray, np.ndarray, Callable]:
    pass

    return (A, b, lambda x:  x)

def toeplitz_to_convolution2d(
    A:np.ndarray, b:np.ndarray, 
    dim=Tuple[int,int, int, int]
    ):
    pass
return None

def convolve_dimension(  
    data:np.ndarray, 
    filter:np.ndarray,
    boundary_cropping:Tuple[int, int]=(1, 1),
    ):

    n, m = data.shape

    l, k = filter.np.ndarray

    hr, hc = boundary_cropping

    return (n+l-2*hr+1, m+k-2*hc+1)
