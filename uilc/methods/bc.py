
import math
from typing import Tuple, Literal,Iterable

import numpy as np
from numpy.polynomial import Chebyshev
from scipy import optimize as op

from uilc import PositionArray
from uilc.methods import esc
from uilc.utils.misc import float_eps, system_param_check, csym_index
from uilc.utils.radiation import lambertian

def _D_f(d, alpha, s):
    return (1+ (0.5*alpha+d)**2)**(-s/2 -1) + (1+ (0.5*alpha-d)**2)**(-s/2-1) - 2*((1+d**2)**(-s/2-1))

def _Di_f(x, s, I0, W, H): # same with (I0/H^2 )*D(x/H, W/H, s)
        return I0/(H)**2 *_D_f(x/H, W/H, s)

def get_de(alpha, s, approx=False): #Done
    app = 0.25*alpha +  math.sqrt(2**(2/(s+2)) -1)/6
    if approx:
        return app
    else:
        r = op.root_scalar(lambda d: _D_f(d, alpha, s), x0 = app, bracket =[0, alpha/2], method = "brentq")
        return r.root
def get_dm(alpha, s, de, approx=False):
    app = math.sqrt(2**(2/(s+2))-1)
    if approx:
        return app
    else:
        r = op.root_scalar( lambda d: _D_f(d, alpha, s) + _D_f(alpha/2, alpha, s), x0 = app, bracket=[0,de], method="brentq")
        return r.root
        
def get_xe(s, W, H, approx=False): #Done
    alpha = W/H
    return H*get_de(alpha, s, approx)
def get_xm(s, W, H, xe, approx=False): #Done
    if xe is None:
        xe = get_xe(s, W, H, approx)
    de = xe/H
    alpha = W/H
    return H * get_dm(alpha, s, de, approx)

def _get_r_points(xarr, dx, s, W, H, xe, xm, threshold=0.7, p_thres= 0.05):
    x_ex = []
    state =0
    exist_p = False
    for i, x in enumerate(xarr):
        if math.fabs(x-0) < float_eps:
            continue
        if x > xe:
            np.delete(xarr,i)
        else:
            x_new = _r_region_x(x,s, W, H, xe, xm, p_thres, exist_p)
            if x<xm:
                exist_p = True
            if x_new is None and state ==0:
                xarr[i] = xe
                state = 1
            elif isinstance(x_new, Iterable):
                x_ex.append(x_new[0])
                x_ex.append(x_new[1])

            else:
                if math.fabs(x_new - x) < dx*threshold and dx+x < W/2:
                    np.delete(xarr,i)
                    x_ex.append(xe)
                else:
                    x_ex.append(x_new)
    x_ex = np.array(x_ex)
    x_ex = np.sort(x_ex)
    return xarr, x_ex

#----------------------------------------------------------------
def _dm_ther(s, W, H, th): # expand permitted area
    alpha= W/H
    th = 1/th if th >1 else th
    
    de = get_de(alpha, s)
    dm = get_dm(alpha, s, de)
    c = (1+th)*_D_f(alpha/2, alpha, s)
    D_f = lambda d: _D_f(d, alpha, s) + c

    r = op.root_scalar(D_f, x0 = dm/2, bracket=[0, de], method="brentq")
    return r.root

def _rq_uniform(s, W, H, th, dif_th=False)->Tuple[int, int]:
    th = 1/th if th>1 else th
    de = get_de(W/H, s)
    dm = get_dm(W/H, s, de)

    if dif_th:
        dm = _dm_ther(s, W, H, th)
        th = 0
    
    dm = (1-th)*dm

    r_em = de/dm
    n =2
    n_even = 0
    n_odd =0
    while(n < 2*(r_em) or n<(r_em)):
        if (n >= r_em and n %2 ==0) and n>=4:
            if n_even == 0:
                n_even = n-2
        n += 1

    n_odd = n-2 if n%2 else n-1 

    return (n_even, n_odd)
def _rq_esc(s, W, H, th, dif_th=False)->Tuple[int, int]:
    th = 1/th if th>1 else th
    de = get_de(W/H, s)
    dm = get_dm(W/H, s, de)

    if dif_th:
        dm = _dm_ther(s, W, H, th)
        th = 0
    n = 2
    dn = esc.coefficient(s, n)[0]
    n_even = n_odd =0
    while(dn <= 2*de/(n-1)):
        if n%2: #odd
            if dm/1.2 <= dn:
                n_odd = n
        else:
            if dm/(0.5+th) <= dn:
                n_even = n
        n +=1
        dn = esc.coefficient(s, n)[0]
    
    return (n_even, n_odd)

def search_rq_fill_numbers(
    s, W, H, 
    method:Literal["uniform", "esc"]="esc", 
    threshold = 0.35, diff_thers=False)->Tuple[int, int]:

    if method == "uniform":
        return _rq_uniform(s, W, H, threshold, diff_thers)
    elif method =="esc":
        return _rq_esc(s, W, H, threshold, diff_thers)
    

def fill_rq(s, W, H, method:Literal["uniform", "esc"]="esc", 
    threshold = 0.35, diff_thers=False, even=False):

    n_even, n_odd = search_rq_fill_numbers(s, W, H, method, threshold, diff_thers)

    #print(n_even, n_odd)
    n = n_even if abs(n_even - n_odd) == 1 else n_odd 
    n = n_even if even else n

    if (n_odd == 0) != (n_even==0):
        n = n_even if n_odd ==0 else n_odd

    if n ==0:
        raise RuntimeError("Failed to find \'n\' value.")

    if method == "uniform":
        xe = get_xe(s, W, H)
        return (2*xe/n) *csym_index(n)
    elif method == "esc":
        return H*esc.array(s, n).get_axis_list()

def get_positive(arr):
    arr = np.sort(arr)
    n = len(arr)
    n_i = int((n-1)/2) if n%2 else int(n/2)
    return arr[n_i:]
def get_full_arr(arr, arr_type=0): # 0 even, 1 odd
    if arr_type: #odd
        arr_pos = arr[1:]
        arr_neg = np.sort(-arr_pos)
        arr =np.concatenate([arr_neg, np.array([0]), arr_pos], axis=0 )
    else:
        arr = np.concatenate([-arr, arr], axis=0)
    return arr

def _r_region_q_x(x, s, W, H, xe, xm):
    alpha = W/H
    de = xe/H
    dm = xm/H
    sgn = math.copysign(1, x)
    d = math.fabs(x/H)

    if d > de or d<0:
        raise ValueError(f"Argument 'x' must be in xm <= x < xe < x <= W/2. {x}")
    
    if math.isclose(d, de, abs_tol = float_eps):
        return de
    
    if math.isclose(d, alpha/2, abs_tol=float_eps) and d<alpha/2:
        return sgn*H*dm
    
    if d < dm: # P region
        raise ValueError(f"x must be larger than {xm}")
    
    if d < alpha/2: # Q,R region
        L = _D_f(d, alpha, s)
        op_f = lambda xc: _D_f(xc, alpha, s) + L
        try:
            sol = op.root_scalar(op_f, bracket=[0, alpha/2], method="brentq")
        except ValueError:
            raise ValueError("x: {}, xe:{}, xm:{}".format(x, H*de, H*dm))
        return sol.root * H * sgn
    
def _get_r_points_q(xarr, dx, s, W, H, xe, xm, thres):
    x_ex = []
    state =0
    for i, x in enumerate(xarr):
        if x > xe:
            np.delete(xarr,i)
        else:
            x_new = _r_region_q_x(x,s, W, H, xe, xm)
            if x_new is None and state ==0:
                xarr[i] = xe
                state = 1
            else:
                if math.fabs(x_new - x) < dx*thres and dx+x < W/2:
                    np.delete(xarr,i)
                    x_ex.append((x+x_new)/2)
                else:
                    x_ex.append(x_new)
    x_ex = np.array(x_ex)
    x_ex = np.sort(x_ex)
    return x_ex
def _get_r_points_p(arr_p, arr_pre, len_arr, d, s, W, H, xe, xm):
    # calculate central radiation.
    arr_pq = PositionArray.from_arrays(get_full_arr(arr_pre, arr_pre.size%2))
    radiation_x = lambda x: lambertian(0, np.sqrt(x), s, H)

    x= np.linspace(-W/2, W/2, int(2*W/d), endpoint=True)
    y = arr_pq.intensity_on(x, radiation_x)
    y_max =y.max()
    boundary_data = y_max*np.ones(y.shape) - y
    
    index_x = np.where(x > xe)
    inter_x = x[index_x:]
    inter_y = boundary_data[index_x:]

    cheby_approx_boundary = Chebyshev.fit(inter_x, inter_y)
    approx_f = Approx2dim()

    srch_domain = [xe, W/2]
    points = []
    while(remained.max() < ):
        
        approx.append(partial(labmertian, t= , s=s, H=H))
        points.append( point )
        srch_domain[1] = point
        pass
    pass


def expansion(arr_pq, s, W, H, thres):
    len_arr = len(arr_pq)
    arr_pq = np.sort(arr_pq)
    d = math.fabs(arr_pq[1]- arr_pq[0])
    pos_arr = get_positive(arr_pq)

    xe = get_xe(s, W, H)
    xm = get_xm(s, W, xe)

    # Sperate to P, Q arrays
    arr_p = arr[arr<=xm]
    arr_q = arr[(arr<=xe) & (arr>xm)]

    arr_rq = _get_r_points_q(arr_p, d, s, W, H, xe, xm)
    arr_pre = np.concatenate([arr_p, arr_q, arr_rq], axis=0)
    arr_rp = _get_r_points_p(arr_p, arr_pre, len_arr, d, s, W, H, xe, xm) # find remained points

    arr_pos = np.concatenate([arr_pre, arr_rp], axis =0)
    arr = get_full_arr(arr_pos, arr_type= len_arr%2)
    return arr




def bc_expansion_2d(arr_pq:PositionArray, s, W:Tuple[float, float], H):
    x_arr = arr_pq.get_axis_list(axis="x")
    y_arr = arr_pq.get_axis_list(axis="y")

    Wx, Wy = W

    x_arr_bc =bc_expansion(x_arr, s, Wx, H)
    y_arr_bc =bc_expansion(y_arr, s, Wy, H)

    return PositionArray.from_arrays(x_arr_bc, y_arr_bc)
