import math
from typing import Tuple, Literal, Callable

import numpy as np
from scipy.optimize import least_squares, nnls

from uilc import PositionArray
from uilc.utils.hypergeo import hyper2F1
from uilc.utils.misc import half_ceil, d2
from uilc.radiation import gaussian, lambertian


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
def propagation_matrix(n:int, W:float, radiation:Callable[[int, int], float])->np.ndarray:
    d= W/n
    return np.fromfunction(lambda i, j: radiation(d*i, d*j), (n,n), dtype=float) 

def _solve_discretized(n, s, W, H, getfd=False):
    F = propagation_matrix(n, W, lambda i,j: lambertian(s, H, d2(i, j)))
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

def solve_nnls(s, W, H, n_nnls = 1, mean=True): #linear, nnls
    d= W/n_nnls
    F = propagation_matrix(n, W, lambda i,j: lambertian(s, H, d2(i, j)))
    delta = nnls(F,np.ones(n_nnls))[0]
    position =np.array([-W/2+d/2+(i-1)*d for i in range(1,n_nnls+1)])
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

# FM method
def _get_freq():
    pass