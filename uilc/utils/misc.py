import math
import numpy as np

float_eps = 1e3*np.finfo(float).eps

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

def extend_signal(sig, time, n=3, period=None, central=True):# not finished
    # Argument check
    if sig.shape != time.shape:
        raise ValueError("Signal and Time dimensions are not same.\n sig:{}, time:{}".format(sig.shape, time.shape))
    if len(sig) > 1:
        raise ValueError("Only 1 dim arrays are permitted. Current:{}".format(sig.shape))
    # Period set
    if period is None:
        period = time[-1]-time[0]
    # Unit signal setting
    time_range = time[-1]-time[0]
    dt = time[1]- time[0]

    if period < time_range:
        if central:
            pass
        else:
            sig = sig[:period]
            time = time[:period]
    else:
        dn = period - time_range
        if central:
            fn = dn/2 if dn%2 ==0 else (dn-1)/2
            bn = dn - fn
            front_sig = np.zeros(shape=(fn,))
            back_sig = np.zeros(shape=(bn,))
            sig = np.concatenate([front_sig +sig +back_sig], axis=0)
        else:
            back_sig = np.zeros(shape=(dn,))
            sig = np.concatenate([sig +back_sig], axis=0)

    sigs = []
    times = []
    if central:
        for i in range(0, n):
            sigs.append(sig + )
    else:
        pass
    

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
