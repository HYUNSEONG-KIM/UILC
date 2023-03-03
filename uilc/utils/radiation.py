import numpy as np
from uilc.utils.misc import float_eps

# Basic functions
def gaussian(x, mu, rho):
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

def lambertian(x, t, s, H):
    k = (x-t)/H
    base = (H**2)*(1 + k**2)**(s/2 +1)
    return np.power(base, -1)
def lambertian2d(r0, r1, s, H):
    if len(r1.shape) >=3 and r1.shape[0] ==2:
        row = r1[0]
        column = r1[1]

        r_v = r0[0] * np.ones(shape=row.shape)
        c_v = r0[1] * np.ones(shape=column.shape)

        k = (np.power(row - r_v, 2) + np.power(column-c_v, 2))/(H**2)
    else:
        k = np.linalg.norm(r0-r1)/H**2
        
    base = (H**2)*(1 + k**2)**(s/2 +1)
    return np.power(base, -1)

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