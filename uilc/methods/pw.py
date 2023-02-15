import math
from typing import Tuple, Literal, Callable

from copy import deepcopy

import numpy as np
from scipy.optimize import least_squares, nnls
from scipy.fft import fft, fftfreq, ifft
from scipy import signal

from uilc import PositionArray
from uilc.utils.hypergeo import get_mpmath_hyper
from uilc.utils.misc import half_ceil, d2, extend_signal
from uilc.radiation import gaussian, lambertian


hyper2F1 = get_mpmath_hyper()

# Math utils
def diff_matrix(n:int):
    return np.array([[1 if i == j else (-1 if i-j == 1 else 0) for j in range(0,n)] for i in range(0,n)])
def inv_diff_matrix(n:int):
    return np.array([[1 if i >= j else (0) for j in range(0,n)] for i in range(0,n)])

def loc_to_diff(x:np.ndarray, n:int):
        fd, = diff_matrix(n)
        d= fd.dot(x)
        w= -2*d[0]
        d0 = d[1:]
        m = math.floor(n/2) if n%2 == 0 else math.floor((n-1)/2)
        return d0[0:m], w, m

def f_residual(d, **kwargs):
        W = kwargs["W"]
        h = kwargs["h"]
        n = kwargs["n"]
        xdata = kwargs["xdata"]
        ydata = kwargs["ydata"]
        m = d.size
        infd = inv_diff_matrix(n)
    
        if n%2 ==0:
            d0 = np.append(d[0:m], np.flip(d[0:m-1]))
            pass
        else:
            d0 = np.append(d,np.flip(d))
    
        d0 = np.insert(d0, 0, -W/2)
        xi = infd.dot(d0)
    
        location = np.array([[[x, 0] for x in xi]])
    
        return (ydata - gaussian(xdata, np.array([[[0]]]), location,h))[0][0]

# Radiation matrix
def position_array(W, n):
    d = W/n
    return np.array([-W/2+d/2+(i-1)*d for i in range(1,n+1)])
def propagation_matrix(n:int, W:float, radiation:Callable[[int, int], float])->np.ndarray:
    d= W/n
    return np.fromfunction(lambda i, j: radiation(d*i, d*j), (n,n), dtype=float) 

def _solve_discretized(n, s, W, H, getfd=False):
    F = propagation_matrix(n, W, lambda i,j: lambertian(s, H, d2(i-j,0)))
    delta = np.linalg.solve(F, np.ones(n))
    if getfd:
        return delta, F
    return delta.min()

def _nmax_app(s:float, W:float, H:float)->int:
    # Get boundary value of positiive solution
    return half_ceil(3/ hyper2F1(1/2, (s+2)/2, 3/2, - (W/(2*H))**2))

def nmax(s:float, W:float, H:float)->int:
    n_app = _nmax_app(s, W, H)
    n =  n_app
    state = 0
    termination = False
    while(not termination):
        minvalue = _solve_discretized(n, s, W, H)
        if minvalue <0:
            if state ==2:
                n -= 1
                termination = True
            else:
                n -=1
                state =1
        else:
            if state ==1:
                termination =True
            else:
                n +=1
                state=2
    return n

def power_weight(s, W, H, dim = 1, set_nmax = True, mean=False): #linear, nnls
    n_max = nmax(s, W, H)
    if set_nmax:
        dim = n_max
    n = dim
    d= W/n
    F = propagation_matrix(n, W, lambda i,j: lambertian(s, H, d2(i-j, 0)))
    
    if n > n_max:
        delta = nnls(F,np.ones(n))[0]
    else:
        delta = np.linalg.solve(F, np.ones(n))

    position =np.array([-W/2+d/2+(i-1)*d for i in range(1,n+1)])
    #Get meaningful points
    if mean:
        therhold = 0.01 * delta.max()
        #therhold = 2* delta.mean()
        position = position[np.argwhere(delta>therhold )[0:,0]]
        delta = delta[delta > therhold ]
    return delta, position, F


#------------------------------------------------------------------------------
# Binarization

# KDE method with gaussian
def nomarlization_lq(arr, xdata, ydata, n, h=False, W=False):
    d, w, m = loc_to_diff(PositionArray.get_axis_list(arr), n)
    infd = inv_diff_matrix(n)
    if h == False:
        h =  w/(1.8*n)
    if W == False:
        W = w
    kwargs = {
                "W": W, 
                "h": h, 
                "n": n, #number of leds in bc solution
                "xdata": xdata, 
                "ydata": ydata
             }
    bc_x0 = d
    sol = least_squares(f_residual, x0 =bc_x0, kwargs=kwargs)
    if sol.status >=1:
        
        if n%2 ==0:
            darr =np.append(sol.x, np.flip(sol.x[0:m-1]))
        else:
            darr = np.append(sol.x, np.flip(sol.x))
        darr = np.insert(darr, 0, -W/2)
        xarr = infd.dot(darr)
        sol_xarr = np.array([[[x, 0] for x in xarr]])
        err  = (sol_xarr[0][n-1][0]-sol_xarr[0][0][0])/2 - sol_xarr[0][n-1][0]
        sol_xarr[0,0:n][0:n,0] = sol_xarr[0,0:n][0:n,0] + err
    else:
        return False
    return sol_xarr


def signal_decomposition(sig, pos, thersholf_freq, return_type=True):
    N = len(pos)
    T = pos[1]-pos[0]

    sig_f = fft(sig)
    x_f = fftfreq(N ,T)

    sig_low = deepcopy(sig_f)
    sig_low[np.where( np.fabs(x_f) > thersholf_freq)] =  0
    sig_high = deepcopy(sig_f)
    sig_high[np.where( np.fabs(x_f) <= thersholf_freq)] =  0

    if return_type:
        return ifft(sig_low), ifft(sig_high)
    else:
        return (sig_low.real, x_f), (sig_high.real, x_f)

def resample_n(sig, time, N, rate=240):
    W = time.max()-time.min()
    sig_ext, time_ext = extend_signal(sig, time, n=N, period =None, central=True)
    sig_resample = signal.resample(sig_ext, N*rate)
    time_resample = position_array(N*W, N*rate)
    return sig_resample, time_resample

def signal_decomposition(sig, pos, thersholf_freq, return_type=True):
    N = len(pos)
    T = pos[1]-pos[0]

    sig_f = fft(sig)
    x_f = fftfreq(N ,T)

    sig_low = deepcopy(sig_f)
    sig_low[np.where( np.fabs(x_f) > thersholf_freq)] =  0
    sig_high = deepcopy(sig_f)
    sig_high[np.where( np.fabs(x_f) <= thersholf_freq)] =  0

    if return_type:
        return ifft(sig_low), ifft(sig_high)
    else:
        return (sig_low.real, x_f), (sig_high.real, x_f)

def get_signal_decomposition(N, s, W, H, ext_n, rate):
    delta, pos, K = power_weight(s, W, H, N, set_nmax=False)
    delta = delta/delta.max()
    sig_ext, pos_ext =resample_n(delta, pos ,ext_n, rate)
    a, b = signal_decomposition(sig_ext, pos_ext, 2*np.pi/W, return_type=False) 
    sig_low, xf1 =a 
    sig_high, xf2 = b 
    return xf1, sig_low, sig_high

#Extend function have a minor centering issue
# Fix tha "uilc.utils.misc.extend_signal"
#def get_peaks(sig, pos, distance=1):
#    index, _ = find_peaks(sig, distance=distance)
#    return np.array([])



