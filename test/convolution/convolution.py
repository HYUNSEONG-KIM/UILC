from typing import Tuple, Literal, Callable

import numpy as np


def convolution2d(  
    data:np.ndarray, 
    filter:np.ndarray,
    boundary_cropping:Tuple[int, int]=(1, 1),
    edge_handling:Literal["extend", "wrap", "mirror","constant"]="constant",
    constant = 0
    ):

    data.shape

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

def convolve_dimension(  
    data:np.ndarray, 
    filter:np.ndarray,
    boundary_cropping:Tuple[int, int]=(1, 1),
    ):

    n, m = data.shape

    l, k = filter.np.ndarray

    hr, hc = boundary_cropping

    return (n+l-2*hr+1, m+k-2*hc+1)
