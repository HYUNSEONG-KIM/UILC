import math
from typing import Tuple, Literal, Callable, Union

from copy import deepcopy

import numpy as np
from scipy.optimize import least_squares, nnls
from scipy.fft import fft, fftfreq, ifft
from scipy import signal

from uilc import PositionArray
from uilc.utils.hypergeo import get_mpmath_hyper
from uilc.utils.mild_math import half_ceil, extend_signal 
from uilc.utils.mild_math import diff_matrix, inv_diff_matrix
from uilc.utils.radiation import gaussian, lambertian


hyper2F1 = get_mpmath_hyper()

# Math utils
def diff_matrix(n:int):
    return np.array([[1 if i == j else (-1 if i-j == 1 else 0) for j in range(0,n)] for i in range(0,n)])
def inv_diff_matrix(n:int):
    return np.array([[1 if i >= j else (0) for j in range(0,n)] for i in range(0,n)])
# Radiation matrix
def position_array(W, n):
    d = W/n
    return np.array([-W/2+d/2+(i-1)*d for i in range(1,n+1)])
def propagation_matrix(n:int, W:float, radiation:Callable[[int, int], float])->np.ndarray:
    d= W/n
    return np.fromfunction(lambda i, j: radiation(d*i, d*j), (n,n), dtype=float) 

def _solve_discretized(n, s, W, H, getfd=False):
    F = propagation_matrix(n, W, lambda i,j: lambertian( i, j, s, H))
    delta = np.linalg.solve(F, np.ones(n))
    if getfd:
        return delta, F
    return delta.min()

def _nmax_app(s:float, W:float, H:float)->int:
    # Get boundary value of positiive solution
    # \frac{3}{{}_2F_1\left(\frac{1}{2}, \frac{s+2}{2}, \frac{3}{2}, - \frac{W^2}{4*H^2}\right)}

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

def power_weight(s, W, H, dim = 1, set_nmax = False, mean=False): #linear, nnls
    n_max = nmax(s, W, H)
    if set_nmax:
        dim = n_max
    n = dim
    d= W/n
    F = propagation_matrix(n, W, lambda i,j: lambertian( i, j, s, H))
    
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

# KDE method with gaussian kernel
def _lq_loc_to_diff(x:np.ndarray, n:int):
        fd, = diff_matrix(n)
        d= fd.dot(x)
        w= -2*d[0]
        d0 = d[1:]
        m = math.floor(n/2) if n%2 == 0 else math.floor((n-1)/2)
        return d0[0:m], w, m

# Gaussian KDE, object function
def _lq_f_residual(d:np.ndarray, **kwargs):
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

def lq_binarization(
    arr:np.ndarray, 
    xdata:np.ndarray, 
    ydata:np.ndarray, 
    n:int, h=False, # Gaussian kernel parameters 
    W=False)-> np.ndarray:

    d, w, m = _lq_loc_to_diff(arr, n)
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
    sol = least_squares(_lq_f_residual, x0 =bc_x0, kwargs=kwargs)
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


def signal_cutting(time, sig, region=[None, None]):
    d_i, d_f = region
    if d_i is not None:
        t_index =  np.argwhere(d_i< time).reshape(-1)
        time = time[t_index]
        sig = sig[t_index]
    if d_f is not None:
        t_index =  np.argwhere(time< d_f).reshape(-1)
        time = time[t_index]
        sig = sig[t_index]
    return time, sig

def signal_decomposition(
    time:np.ndarray, 
    sig:np.ndarray, 
    thersholf_freq:float, 
    return_type=True)-> Tuple[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]]:

    N = len(time)
    T = time[1]-time[0]

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

def resample_n(time:np.ndarray, sig:np.ndarray, N:int, rate=240):
    W = time.max()-time.min()
    sig_ext, time_ext = extend_signal(sig, time, n=N, period =None)
    #sig_resample, time_resample  = signal.resample(sig_ext, N*rate, t= time_ext)
    sig_resample  = signal.resample(sig_ext, N*rate)
    time_resample = position_array(N*W, N*rate) - (0 if N%2 else W/2)

    try:
        time_resample -= detecting_shifted_time(time_resample, sig_resample, [-W/2, W/2])
    except IndexError:
        try:
            time_resample -= detecting_shifted_time(time_resample, sig_resample, [-W, W])
        except:
            pass
    return sig_resample, time_resample

def detecting_shifted_time(
    time:np.ndarray, sig:np.ndarray, 
    domain_range=[None,None])->float: # Shifted even function restoring.
    # 1. extream detection
    # 2. ordering by the distance from 0 point
    time, sig = signal_cutting(time, sig, region = domain_range)
    d_i, d_f = domain_range
    if d_i is not None and d_f is not None:
        center = (d_f + d_i)/2
    else:
        center = 0
    index, _ = signal.find_peaks(sig)

    x_peaks = time[index]
    peaks = sig[index]

    index_neg = np.argwhere(x_peaks<center).reshape(-1) #index negative
    index_pos = np.argwhere(x_peaks>center).reshape(-1) #index positiive

    i_n1 = index_neg[-1]
    i_n2 = index_neg[-2]
    i_p0 = index_pos[0] 
    i_p1 = index_pos[1]

    p_n2 = peaks[i_n2]
    p_n1 = peaks[i_n1]
    p_p0 = peaks[i_p0]
    p_p1 = peaks[i_p1]

    d1 = math.fabs(p_p0 - p_n1)
    d2 = math.fabs(p_p1 - p_n1)
    d3 = math.fabs(p_p0 - p_n2)
    
    b_pos = d2<d1
    b_neg = d3<d1
    if len(index_neg) < len(index_pos):
        b_neg = False
    elif len(index_neg) > len(index_pos):
        b_pos = False

    if b_pos == b_neg:
        i0 = i_p0
        i1 = i_n1
    elif b_pos:
        i1 = i0 = i_p0
    else:
        i1 = i0 = i_n1
        
    center_shifted = (x_peaks[i0] + x_peaks[i1])/2

    return center_shifted

def get_signal_decomposition(N, s, W, H, ext_n, rate):
    delta, pos, K = power_weight(s, W, H, N, set_nmax=False)
    delta = delta/delta.max()
    sig_ext, pos_ext =resample_n(pos, delta, ext_n, rate)
    
    a, b = signal_decomposition(
        pos_ext, sig_ext,
        2*np.pi/W, return_type=False) 
    sig_low, xf1 =a 
    sig_high, xf2 = b 
    return xf1, sig_low, sig_high



def fourier_binarization(
    time:np.ndarray, weight:np.ndarray, 
    ext_n:int, rate:int,
    height:Tuple[Union[None, float], Union[None, float]]=[None, None],
    threshold:Tuple[Union[None, float], Union[None, float]]=[None, None],
    distance:Tuple[Union[None, int], Union[None, int]]=[None, None]
    ):
    weight = weight/weight.max()
    W = (time.max()-time.min()) + (time[1]-time[0])

    sig_ext, pos_ext = resample_n(time, weight, ext_n, rate)
    sig_low_de, sig_high_de = signal_decomposition(sig_ext, pos_ext, 2*np.pi/W, return_type=True)
    
    #time_low, sig_low = signal_cutting(pos_ext, sig_low_de.real, region = [-W/2 , W/2])
    #time_high, sig_high = signal_cutting(pos_ext, sig_high_de.real, region = [-W/2 , W/2])
    #
    #return time_low, sig_low, time_high, sig_high
    #peaks_sig_low, _ = signal.find_peaks(sig_low, height=height[0], threshold=threshold[0], distance= distance[0])
    #peaks_sig_high, _ = signal.find_peaks(sig_high, height=height[1], threshold=threshold[1], distance= distance[1])
#
    #return time_low[peaks_sig_low], time_high[peaks_sig_high]




    


