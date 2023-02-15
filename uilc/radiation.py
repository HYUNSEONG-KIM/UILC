import numpy as np
from uilc.utils.misc import float_eps

def lambertian_1d(xarr, sources, s, W, H):
    result = np.zeros(xarr.shape)
    for s in sources:
        result += H**s/((s - xarr)**2 + H**2)**(s/2+1)
    return result
def lambertian(s, h, d, inv=True):
    r = (h**2 + d)**-1 if inv else 1
    return r/(1 + d/(h**2))**(s/2)
def gaussian(s, h, d, inv = False):
    r =  (1/(h**2 + d)) if inv else 1.
    return r*np.exp(- s*(np.sqrt(d)/h)**2)
def gauss_tan(s, h, d, inv=False):
    r =  (1/(h**2 + d)) if inv else 1.
    return r*np.exp(- s*(np.tan(np.abs(np.sqrt(d))/h)**2))
def dirac(d, ep=1000*float_eps):
    return ((np.sqrt(np.pi)*ep))**(-1) * np.exp(- d/(ep)**2)