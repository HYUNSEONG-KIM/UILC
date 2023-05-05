from typing import Tuple, Union, Iterable
from numbers import Number

import numpy as np
from numpy.typing import NDArray
from uilc.utils.misc import float_eps

# Basic functions
def gaussian_kernel(x, mu, rho):
    return (1/(rho * np.sqrt(2 * np.pi)))*np.exp(- 0.5*((x-mu)/rho)**2)
def dirac(x, ep=1E3*float_eps):
    # Gaussian dirac sequence approximation
    return ((np.sqrt(np.pi)*ep))**(-1) * np.exp(- (x/ep)**2)
def kronecker_delta(i, j):
    return 1 if i==j else 0

# On plane
# x: source plane
# t: above plane
def inverse_law(x, t, H):
    return (H**2 + (x-t)**2)**-1

def lambertian(x:Union[Number, NDArray], t:float, s:float, H:float)-> float:
    """Lambertian radation function of 1 dimension. Calculating `s` power radiation of `x` sources to `t` point on `H` distance plane. 

    Args:
        x (Union[Number, NDArray]): Point(s) of source(s).
        t (float): Point of target plane.
        s (float): Lambertian radiation power parameter.
        H (float): Distance from sources and target plane.

    Returns:
        float: Radiation value of the given `(t, H)` point. 
    """
    k = (x-t)/H
    base = (H**2)*(1 + k**2)**(s/2 +1)
    return (1/base).sum() if isinstance(x, np.ndarray) else (1/base)

def lambertian2d(sources:NDArray, target_plane:Tuple[NDArray, NDArray], s:float, H:float) -> NDArray:
    """Lambertian radation of power `s` calculation for the given `sources` distribtuon on `target_plane`.

    Args:
        sources (NDArray): [[x1, x2, ...], [y1, y2, ...]]. Meshgrid [X, Y] will automatically transformed to presented form.
        target_plane (Tuple[NDArray, NDArray]): X, Y meshgrid of the target plane
        s (float): Lambertian radiation power parameter.
        H (float): Distance from sources and target plane.

    Returns:
        NDArray: Radiation value on the target plane, which have a same dimension with X, Y array of `target_plane`
    """
    if (not isinstance(sources, np.ndarray)) and isinstance(sources, Iterable):
        sources = np.array(sources)
    elif isinstance(sources, Number):
        sources = np.array([sources, 0])
    
    if not (isinstance(sources, np.ndarray) or isinstance(sources, Number)): 
        raise TypeError("Sources must be a ndarray or at least be a number.")

    if len(sources.shape) ==3 and sources.shape[0] == 2: # Meshgrid
        X = sources[0].reshape(-1)
        Y = sources[1].reshape(-1)
        sources = np.vstack([X, Y])
    if sources.shape[-1] == 2 and len(sources.shape) >1: # Set of (x,y) points
        sources = sources.reshape((-1, 2))
    
    x, y = sources
    X_t, Y_t = target_plane

    if isinstance(x, Iterable):
        x_d = np.stack(X_t - xi for xi in x)
        y_d = np.stack(Y_t - yi for yi in y)
        D = x_d**2 + y_d**2
    else:
        D = (X_t-x)**2 + (Y_t-y)**2
    
    D = D/ H**2
    base = H**2 * (1+ D)**(s/2+1)
    cal_result = np.power(base,-1)

    if len(cal_result.shape) >2 :
        return cal_result.sum(axis=0)
    else:
        return cal_result


def gauss_tan_radi(x, t, s, h):
    r =  (1/(h**2 + (x-t)**2))
    return r*np.exp(- s*(np.tan((x-t)/h)**2))

def gauss_radi(x, t, s, h, inv = False):
    r =  (1/(h**2 + (x-t)**2)) if inv else 1.
    return r*np.exp(- s*((x-t)/h)**2)

def lambertian_convolve_1d(
    xarr:np.ndarray, 
    sources:np.ndarray, 
    s:float, 
    H:float):
    result = np.zeros(xarr.shape)
    for s in sources:
        result += H**s/((s - xarr)**2 + H**2)**(s/2+1)
    return result