from typing import Tuple, Callable
#------------------------
import math

import numpy as np
from scipy.fft import fft, fftfreq, ifft
from scipy import signal
#-----------------------

def r2(x:float,y:float)->float:
    return x**2 + y**2

# Integer - Number utils---------------------------------------------
def half_ceil(x:float)->int:
        rho = x - math.floor(x)
        if rho >= 0.5:
            result = math.ceil(x)
        else:
            result = math.floor(x)
        return result
def near_integers(x:float)->int:
        low = math.floor(x)
        return low, low+1 
# Matrix -vector routines---------------------------------------------
def arithmetic_sequence(N:int)->np.ndarray: # = np.arange(N)+1
    return np.array([i+1 for i in range(0,N)])
def center_sym_index(N:int)->np.ndarray:
    return np.array([(i-(N-1)/2) for i in range(0, N)])
def array2meshgrid(
    x_range:np.ndarray, 
    y_range:np.ndarray, 
    dim:int)->np.ndarray: #Convert 2 1 dim arrays to meshgrid
    xline = np.linspace(x_range[0], x_range[1], dim[0])
    yline = np.linspace(y_range[0], y_range[1], dim[1])
    return np.meshgrid(xline, yline)
def region_mesh(
    boundary:Tuple[float, float, float, float], 
    dim:int,
    scale:float=1)->Tuple[np.ndarray, Tuple[float, float, float, float]]:
        scale_x, scale_y = scale
        xi, xf, yi, yf = boundary
        x_r = (scale_x * xi, scale_x * xf)
        y_r = (scale_y * yi, scale_y * yf)
        return array2meshgrid(x_r, y_r, dim), [*x_r, *y_r] 

# Matrix
def diff_matrix(n:int)->np.ndarray:
    return np.array([[1 if i == j else (-1 if i-j == 1 else 0) for j in range(0,n)] for i in range(0,n)])
def inv_diff_matrix(n:int)->np.ndarray:
    return np.array([[1 if i >= j else (0) for j in range(0,n)] for i in range(0,n)])

def position_array(W, n):
    d = W/n
    return np.array([-W/2+d/2+(i-1)*d for i in range(1,n+1)])

# Radiation matrix
def propagation_matrix(
    n:int, W:float, 
    radiation:Callable[[int, int], float])->np.ndarray:
    d= W/n
    return np.fromfunction(lambda i, j: radiation(d*i, d*j), (n,n), dtype=float) 

# Data 
def data_ceiling(data:np.ndarray, n:int)->np.ndarray:
    if n <= 1:
        raise ValueError("n must be greater than 1 and integer.")
    d_i = data.min()
    dn = (data.max()-data.min())/n
    dn_add = (data.max()-data.min())/(n-1)
    for i in range(0, n):
        di1 = i*dn + d_i
        di2 = di1+dn
        data = np.where((data >= di1) & (data < di2 ), ((i))*dn_add + d_i, data)
    return data

# Signal routines
def extend_signal(
    sig:np.ndarray, 
    time:np.ndarray, 
    n:int=3, 
    period:float=None, 
    central:bool=True, 
    edge_matching:bool = False)-> Tuple[np.ndarray, np.ndarray]:# not finished
    # Argument check
    if sig.shape != time.shape:
        raise ValueError("Signal and Time dimensions are not same.\n sig:{}, time:{}".format(sig.shape, time.shape))
    if len(sig.shape) > 1:
        raise ValueError("Only 1 dim arrays are permitted. Current:{}".format(sig.shape))
    # Period set
    if period is None:
        period = time[-1]-time[0]
    # Unit signal setting
    time_range = time[-1]-time[0]
    dt = time[1]- time[0]
    period_n = half_ceil(period/dt)

    if period < time_range:
        if central:
            resi = sig.shape[0] - period_n
            k = int((resi+1)/2) if resi %2 ==1 else int(resi/2)

            sig = sig[k: k+period_n]
            time = time[k: k+period_n]
        else:
            sig = sig[:period_n]
            time = time[:period_n]
    else:
        dn = period_n - time.shape[0]
        if dn <0 :
            dn = 0
        if central:
            fn = int(dn/2) if dn%2 ==0 else int((dn-1)/2)
            bn = dn - fn

            front_sig = np.zeros(shape=(fn,))
            back_sig = np.zeros(shape=(bn,))
            sig = np.concatenate([front_sig, sig, back_sig], axis=0)
            
            front_t = np.array([ time.max() + dt*(i +1) for i in range(0, fn)])
            back_t = np.array([ time.min() - dt*(bn - i -1) for i in range(0, bn)])
            time = np.concatenate([front_t, time, back_t], axis=0)
        else:
            back_sig = np.zeros(shape=(dn,))
            sig = np.concatenate([sig, back_sig], axis=0)
            back_time = [dt*(i+1) + time.max() for i in range(0, dn) ]
            time = np.concatenate([time, back_time], axis=0)

    # Extend the given signal to n-th period signal.
    time = np.sort(time)
    sigs = []
    times = []

    for i in range(0, n):
        n_i = i+1 - n/2
        time_i = time + n_i * period
        times.append(time_i)
        sig_i = sig[::-1] if edge_matching and i%2 else sig
        sigs.append(sig_i)

    
    return np.concatenate(sigs, axis=0), np.sort(np.concatenate(times, axis=0))
    
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
    return_type=True
    )-> Tuple[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]]:
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

def detecting_shifted_time_max(time:np.ndarray, sig:np.ndarray):
    sig_max = sig.max()

    index_peaks, _ = signal.find_peaks(sig, height = 0.95*sig_max)

    if len(index_peaks) <2:
        return None
    else:
        i1 = index_peaks[0]
        i2 = index_peaks[1]

        return (time[i1] + time[i2])/2

def resample_n(time:np.ndarray, sig:np.ndarray, N:int, rate=240):
    W = time.max()-time.min()
    sig_ext, time_ext = extend_signal(sig, time, n=N, period =None)
    sig_resample  = signal.resample(sig_ext, N*rate)
    time_resample = position_array(N*W, N*rate) - (0 if N%2 else W/2)

    time_resample, sig_resample = signal_cutting(time_resample, sig_resample, [-W, W])

    dt = detecting_shifted_time_max(time_resample, sig_resample)
    if dt is None:
        dt = 0

    time_resample = time_resample - dt
    return sig_resample, time_resample