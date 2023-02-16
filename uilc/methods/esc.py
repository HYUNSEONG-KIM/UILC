import math
from typing import Tuple, Literal, Union

import numpy as np
from scipy.optimize import root_scalar, minimize_scalar

from uilc import PositionArray

from uilc.utils.misc import float_eps # constants
from uilc.utils.mild_math import r2, half_ceil, center_sym_index
from uilc.utils import bresenham # module

# Expanded Sparrow Criterion
# Implementation and additional utils
# Reference: I. Moreno et al, 2006 , Designing light-emitting diode arrays for uniform near-field irradiance, Appl. Opt., vol 45, 10, p 2265, 4 2006, doi: 10.1364/AO.45.002265.

def _parameter_check(s:float, N:int, M:int):
    if s <1.0:
        raise ValueError(f"\'s\' must be 1 or more real number. \n s:{s}, {type(s)}")
    if N * M == 0 or N*M < float_eps:
        raise ValueError(f"\'N\' and \'M\' must be 1 or more integers. \n N:{N}, {type(N)}\n M:{M}, {type(M)}")

def _op_func_linear(D:float, s:float, N:int)->float: # _linear
    y =0.0
    for i in range(1,N+1):
        y += (1-(s+3)*(N+1-2*i)**2 * (D**2)/4)*((N+1-2*i)**2 * (D**2)/4 +1)**(-(s+6)/2)
    return y

def _op_func_rectangularangular(D:float, s:float, N:int, M:int)->float: # _rectangularangular
    y =0.0
    for i in range(1,N+1):
        for j in range(1, M+1):
            y += (((N+1-2*i)**2 + (M+1-2*j)**2)*(D**2/4)+1)**(-(s/2+3.0)) * (1-((s+3)*(N+1-2*i)**2 -(M+1-2*j)**2)*(D**2)/4)
    return y
    
def _coef_linear( s:float, N:int, approx:bool = False)->Tuple[float, float]: #_coefficient_linear
    if N == 2:
        cof = math.sqrt(4/(s+3))
    elif approx and (N >4 and s >30):
        cof = math.sqrt(3.2773/(s+4.2539))
    elif N%2 == 0 :
        sol = root_scalar(lambda D: _op_func_linear(D, s, N), bracket=[0,1], method = "brentq")
        cof = sol.root
    else:
        res = minimize_scalar(lambda D: _op_func_linear(D, s, N), bounds=(0,1), method = "bounded")
        cof = res.x
    return cof
    
def _coef_rectangular( s:float, N:int, M:int, approx:bool=False)->Tuple[float, float]:
    if M > N:
        N, M = M, N
    if N==2 and N == M:
        cof= math.sqrt(4/(s+2))
    if approx == True and (N > 4 and M > 4 and s>30):
        cof= math.sqrt(1.2125/(s-3.349))
    else:
        try:
            sol = root_scalar(lambda D: _op_func_rectangularangular(D, s, N, M), bracket=[0,1],method="brentq")
            if sol.converged == False:
                res = minimize_scalar(lambda D: _op_func_rectangularangular(D, s, N, M), bounds=(0,1), method = "bounded")
                cof= res.x
            else:
                cof= sol.root
        except:
            res = minimize_scalar(lambda D: _op_func_rectangularangular(D, s, N, M), bounds=(0,1), method = "bounded")
            cof= res.x
    return cof

def coefficient(
    s:float,
    N:int, M:int=1,
    shape:Literal["L", "R"]="L",
    approx:bool=False
)-> Tuple[float, float]:
    _parameter_check(s, N, M)
    N = half_ceil(N)
    M = half_ceil(M)
    if shape == "L":
        cof_x = _coef_linear(s, N, approx)
        cof_y = None if M == 1 else _coef_linear(s, M, approx)
    elif shape == "R":
        cof_x = cof_y = _coef_rectangular(s, N, M, approx)
    else:
        raise ValueError("\"shape\" argument must be \"L\" or \"R\" current value is {}".format(shape))
    return cof_x, cof_y

#-------------------------------------------------------------------------------
def area(s:float, nx:int, ny:int=1, shape:Literal["L", "R"]="L"):
    _parameter_check(s, nx, ny)
    dx, dy = coefficient(s, nx, ny, shape=shape)
    if shape == "L":
        return (nx-1)*dx, 0
    elif shape == "R" and ny >=1:
        return (nx -1)*dx, (ny-1)*dy
    else:
        raise ValueError(f"shape:{shape}, ny:{ny}")

# Linear uniform distribution search --------------------------------------
def _linear_nmax( 
        s:float, W:float, H:float, 
        threshold:float=0.3, permit_exceed=True
        )->Union[int, None]:
    xlim =W/H
    n_2 = 2
    n_3 = 3

    length_2 = area(s, n_2)[0]
    length_3 = area(s, n_3)[0]
    
    if length_2 > xlim:
        n = n_2
    elif length_2 < xlim and length_3 > xlim:
        n = n_3
    else:
        n = n_3
        length = length_3
        while (length < xlim):
            n += 1
            length = area(s, n)[0]

    # Choose a closest number to boundary
    if n > 3:
        ns = [n-2, n-1, n, n+1]
    else:
        if n == 3: #check 2, 3
            ns = [2, 3, 4, 5]
        else:
            ns = [2, 3, 4]
    ds = [coefficient(s, n)[0] for n in ns]
    areas = [ area(s, n)[0] for n in ns]
    residual =[(xlim-area)/2 for area in areas]

    # Check exceeding permit
    if not permit_exceed:
        ns = [n for n, resi in zip(ns, residual) if resi > 0]
        ds = [d for d, resi in zip(ds, residual) if resi > 0]
        areas = [area for area, resi in zip(areas, residual) if resi > 0]
        residual =[resi for resi in residual if resi>0]

    # threshold check
    thresholds = [ math.fabs(resi) < threshold*d for d, resi in zip(ds, residual)]
    residual = [math.fabs(resi) for resi in residual]
    residual_ther =[resi for resi, thers in zip(residual, thresholds) if thers]
    residual_not_ther =[resi for resi, thers in zip(residual, thresholds) if not thers]
    n_ther = [ n for n, thers in zip(ns, thresholds) if thers ]
    n_not_ther = [ n for n, thers in zip(ns, thresholds) if not thers ]
    
    if len(n_ther) != 0:
        N = n_ther[np.argmin(residual_ther)]
        ther = True
    elif len(n_not_ther) != 0:
        N = n_not_ther[np.argmin(residual_not_ther)]
        ther = False
    else:
        raise RuntimeError("Impossible")
    return (N, N, 1, ther)

# Rectangular uniform distribution search --------------------------------------

# status
a = 1 # both larger than hor and ver values of the area
b = 0 # larger, smaller or smaller, larger
c = -1 # both smaller than the hor and ver values of the area        
def _rect_point_clasification(p, s, W):
    Wx, Wy = W
    nx, ny = p
    area_ = area(s, nx, ny, shape="R")
    x_b, y_b = area_[0] > Wx, area_[1] > Wy
    if x_b == y_b:
        if x_b and y_b: # > , >
            status = -1
        else: # <, <
            status = 1
    else:
        status = 0
    return status
    
def _rect_point_residual(p, s, W):
    Wx, Wy = W
    nx, ny = p
    area_ = area(s, nx, ny, shape="R")
    return (Wx - area_[0], Wy - area_[1])
    
def _rect_threshold(points, s, W, threshold):
    if points is None:
        return [], []
    residual= [ _rect_point_residual(p, s, W) for p in points ]
    ds = [ coefficient(s, *p, shape="R")[0] for p in points ]

    point_stricted =[ p for p, resi, d in zip(points, residual, ds) if math.fabs(resi[0]) < d*threshold and math.fabs(resi[1]) < d*threshold]
    point_unstricted =[ p for p, resi, d in zip(points, residual, ds) if not (math.fabs(resi[0]) < d*threshold and math.fabs(resi[1]) < d*threshold)]

    return point_stricted, point_unstricted

def _rect_exceed_and_threshold( points, s, W, threshold, permit_exceed):
    if not permit_exceed and points is not None:
        residuals=[_rect_point_residual(p, s, W) for p in points]
        points= [p for p, resi in zip(points, residuals) if resi[0]>0 and resi[1]>0]
    return _rect_threshold(points, s, W, threshold)

def _rect_get_point_range(
    initiation:Tuple[Tuple[float, float], int, float], 
    s:float, W:Tuple[float, float], iter_max = 300):
    # Search
    # # 1. Get adjacent points near 'y = mx' line for given 'x'.
    # # 2. Calculate measure of (N_i, M_i) points at 1 -> a, b, c 
    # # 3. Find a- > c brackets near 'y = mx'.
    # # 3. If exceeding event is occured, pick current and previous point.
    Wx, Wy = W
    m = Wy/Wx
    point_range = [None, None]
    p, status, residual = initiation
    dir = True if status != c else False
    p_range = (0, 1) if dir else (1, 0)

    iteration_max = iter_max
    iter_i = 0
    while(iter_i < iteration_max):
        iter_i += 1
        # Get new points and estimate them
        p_current= bresenham.points(p, line_param=(m,0), p_range=p_range)[0][0]
        status_current = _rect_point_clasification(p_current, s, [Wx, Wy])
        residual_current = _rect_point_residual(p_current, s, [Wx, Wy])
        # Break, or manipulate direction
        if status_current != status:
            if status_current == b:
                if status == a:
                    point_range[0] = p
                else:
                    point_range[1] = p
            else: # direction change
                if status == -status_current:
                    pi, pf = (p, p_current) if p_current[0] > p[0] else (p_current, p) 
                    point_range[0] = pi
                    point_range[1] = pf
                    break
                elif status == b:
                    if status_current == a:
                        point_range[0] = p_current
                    else:
                        point_range[1] = p_current
                
                dir = not dir
                p_range = (0, 1) if dir else (1, 0)
                continue
        else: # Check that residual is decreasing or increasing.
            if status != b: # a and c states, same if "status_current != b"
                dfx, dfy = residual_current[0] - residual[0], residual_current[1] - residual[1]
                if dfx*dfy >0 and status * dfx >0: 
                    # ((fx'>0 and fy'>0) or (fx'<0 and fy'<0)) 
                    # and ((a:1 and f'>0) or (c:-1 and f'<0))
                    dir = not dir
                    p_range = (0, 1) if dir else (1, 0)
                    continue
        if point_range[0] is not None and point_range[1] is not None:
            break

        p = p_current
        status = status_current
        residual = residual_current
            
    if iter_i >= iteration_max and point_range[0] is None:
        raise RuntimeError(f"Cannot find point in {iteration_max} iteration.")
    return point_range

def _rectangular_nmax(
    s:float, 
    W:Tuple[float, float], H:float, 
    threshold:float=0.3, permit_exceed=True,
    iter_max = 300
    ):
    # This function search point range 
    # using Bresenham line algorithm, "bresenham" module.

    Wx, Wy = W
    m = Wy/Wx
    d = 0
    Wx, Wy = Wx/H, Wy/H
        

    # initial point setting
    if s>30:
        approx_d = coefficient(s, N=5, M=5, approx=True)
        p_i = [half_ceil(Wx/approx_d), half_ceil(Wy/approx_d)]
    else:
        # Using line nmax
        n_x = _linear_nmax(s, Wx, H, permit_exceed)[1]
        p_i = [n_x, half_ceil(m * n_x)] 
    
    if p_i[1] <2:
        p_i[0] = half_ceil(2/m)
        p_i[1] = 2
    
    residual_i = _rect_point_residual(p_i, s, [Wx, Wy])
    status_i = _rect_point_clasification(p_i, s, [Wx, Wy])
    
    
    point_range = _rect_get_point_range((p_i, status_i, residual_i),s, W=[Wx, Wy], iter_max = iter_max)

    # Get estimation range--------------------------------------------------
    pi, pf = point_range[0], point_range[1]
    points_main, points_sub = bresenham.points(
        pi, line_param = [m, d], 
        p_range=(2, pf[0]-pi[0]+2), allow_cross_p = True)

    # Estimation - in progress
    main_strict, main_unstrict = _rect_exceed_and_threshold(points_main, s, [Wx, Wy], threshold, permit_exceed)
    sub_strict, sub_unstrict = _rect_exceed_and_threshold(points_sub, s, [Wx, Wy], threshold, permit_exceed)
    #--------------------------------------------------------------------------
    threshold_list = main_strict + sub_strict
    unthreshold_list = main_unstrict + sub_unstrict
        
    final_list = threshold_list if len(threshold_list) !=0 else unthreshold_list

    if final_list is None or len(final_list) == 0:
        print("Stricted:-------------------------")
        print(threshold_list)
        print("Unstricted:-------------------------")
        print(unthreshold_list)

        raise RuntimeError("Cannot find appropriate points.")
    else:
        norms = [ r2(*_rect_point_residual(p, s, [Wx, Wy])) for p in final_list]
        nx, ny = final_list[np.argmin(norms)]
        return nx*ny, nx, ny, len(threshold_list) !=0
    
def nmax_for_region(
    s:float, 
    W:Union[float, Tuple[float, float]], H:float, 
    shape:Literal["L", "R"] = "L", 
    threshold:float = 0.3,
    permit_exceed:bool = True
    ) -> Tuple[int, int, int, bool]:

        if hasattr(W, "__iter__"):
            if len(W) >=2:
                W = W[:2]
                Wx, Wy = W
            else:
                Wx =W[0]
                Wy = None
        else:
            Wx = W
            Wy = None

        if shape == "L":
            n1, m1, l1, ther_1 = _linear_nmax(s, Wx, H, threshold, permit_exceed)
            N = (n1, m1, l1, ther_1)
            if Wy is not None:
                n2, m2, l2, ther_2 = _linear_nmax(s, Wy, H, threshold, permit_exceed)
                N = (m1 * m2, m1, m2, ther_1 and ther_2)
        elif shape =="R":
            if Wy is None:
                Wy = Wx
            switch = False
            if Wy > Wx:
                Wx, Wy = Wy, Wx
                switch = True
            nm, n, m, ther = _rectangular_nmax(s, (Wx, Wy), H, threshold, permit_exceed)
            N = (nm, m, n , ther) if switch else (nm, n, m, ther)
        else:
            raise ValueError(f"\"shape\" must be \"L\" or \"R\" current:{shape}")
        return N
    
def array(s:float, N:int, M=1, shape="L", approx=False):
    _parameter_check(s, N, M)

    dx, dy = coefficient(s, N, M, shape, approx)
    if dy is None:
        dy = 1
        
    xarr = dx*center_sym_index(N)
    yarr = dy*center_sym_index(M)
    return PositionArray.from_arrays(xarr, yarr)
      