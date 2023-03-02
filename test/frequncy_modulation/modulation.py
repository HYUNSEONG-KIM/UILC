from scipy.interpolate import interp1d
import numpy as np
from scipy.signal import cwt
from scipy.signal import ricker as wavelet_ricker
from scipy.signal import argrelextrema


#FM binarization

def get_period_ds(position, y):
    index =argrelextrema(y, np.greater)
    period = 0
    for i in index:
        if i+2 < len(y):
            period += (position[i+2]- position[i])
    period = period / len(y)
    return period

def get_period_nnls(position, y):
    # Wavelet transform
    # 
    # Get Extrem value from transformed data
    # 
    # Calculate Period
    pass
def get_base_f_ds(position, y):
    return 1/get_period_ds(position, y)

def get_base_f_nnls(position, y):
    return 1/get_period_nnls(position, y)



if __name__ == "__main__":
    import sys, os
    sys.path.insert(1, os.getcwd())
    sys.path.append("..")

    import src.uilc as uilc

    cm = 1E-2
    s = 30
    W = 9*cm
    H = 3*cm
    ds_rho, ds_position, ds_F = uilc.disop.solve_linear(s, W, H)
    f_b = get_base_f_ds(ds_position, ds_rho)
