import math
from typing import Tuple
from numbers import Number

import numpy as np

float_eps = 1e3*np.finfo(float).eps

def system_param_check(s:float, W:Tuple[float, float], H:float):
    if not isinstance(s, Number):
        raise TypeError(f"s: It must be a numberic type. \nCurrent: {s}, {type(s)}")
    if s <1 :
        raise ValueError(f"s: It must be 1 or more real number.\nCurrent:{s}")
    
    if not hasattr(W, "__iter__"):
        if isinstance(W, Number):
            W = [W]
        else:
            raise TypeError(f"W: It must be a numeric tuple of length 2 or a  numeric variable.\nCurrent: {W}, {type(W)}")
    if len(W) >2:
        raise ValueError(f"W: It must be a numeric tuple of length 2 or a numeric variable.\nCurrent: {W}, {type(W)}")
    if not (isinstance(W[0], Number)) or not (isinstance(W[1], Number)):
        raise ValueError(f"W: It must be a numeric tuple of length 2 or a numeric variable.\nCurrent: {W}, {type(W[0])}, {type(W[1])}")
    if W[0] <0 or W[1] <0:
        raise ValueError("W: It must be a tuple of positive real numbers.\nCurrent: {W}")


    if not isinstance(H, Number):
         TypeError(f"W: It must be a numeric variable.\nCurrent: {H}, {type(H)}")
    if H < 0:
        raise ValueError(f"H: It must be positive real number. \nCurrent: {H}")
    
    return 0

def d2(x,y):
    return x**2 + y**2

def plane_meshgrid(x_range, y_range, dim):
    xline = np.linspace(x_range[0], x_range[1], dim[0])
    yline = np.linspace(y_range[0], y_range[1], dim[1])
    return np.meshgrid(xline, yline)

def plot_xy_mesh(points, scale, dim):
        scale_x, scale_y = scale
        xi, xf, yi, yf = points
        x_r = (scale_x * xi, scale_x * xf)
        y_r = (scale_y * yi, scale_y * yf)
        return plane_meshgrid(x_r, y_r, dim), [*x_r, *y_r] 

def csym_index(N):
    return np.array([(i-(N-1)/2) for i in range(0, N)])

def half_ceil(x):
        rho = x - math.floor(x)
        if rho >= 0.5:
            result = math.ceil(x)
        else:
            result = math.floor(x)
        return result
def near_integers(x):
        low = math.floor(x)
        return low, low+1 
        
def data_ceiling(data, n):
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
    

def rectangle_line_points(dx, dy, Wx=None, Wy=None, wx=1, wy=1):
    Wx = 0 if Wx is None else Wx
    Wy = 0 if Wy is None else Wy 
    return [
        [[wx*(Wx-dx)/2, wx*(Wx+dx)/2],[wy*(Wy-dy)/2, wy*(Wy-dy)/2]],
        [[wx*(Wx-dx)/2, wx*(Wx+dx)/2],[wy*(Wy+dy)/2, wy*(Wy+dy)/2]],
        [[wx*(Wx-dx)/2, wx*(Wx-dx)/2],[wy*(Wy-dy)/2, wy*(Wy+dy)/2]],
        [[wx*(Wx+dx)/2, wx*(Wx+dx)/2],[wy*(Wy-dy)/2, wy*(Wy+dy)/2]]
    ]

def print_mesh_point(arr):
    for line in arr:
        for element in line:
            print(f"({element[0]:.2}, {element[1]:.2})", end = "")
        print(":\n")
