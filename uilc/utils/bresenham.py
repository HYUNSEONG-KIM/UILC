
# Dot searching algrithm using Bresenham line algorithm
# refer:
#   Bresenham line algorithm
#       - https://www.cs.helsinki.fi/group/goa/mallinnus/lines/bresenh.html
#       - https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm, 10, Feb, 2023

import math
from typing import Tuple, Union

import numpy as np

from uilc.utils.misc import float_eps

def line(x, m, b):
    return m * x + b
def plane_deter(point, m, d, dx = 1):
    a = m* dx
    b = - dx
    c = dx* d
    x, y = point
    result = (a*x + b*y + c)
    if math.fabs(result) < float_eps:
        place = 0
    elif result >0:
        place = 1
    else:
        place = -1
    return place
def f_generator(m, d):
    return lambda x, y: plane_deter((x,y), m , d)

# p -> p1, p2 ,p3: p1(diagonal), p2(prime), p3(sub)
def next_points(point, m, dir = True):
    xi, yi = point
    dx = 1 if dir else -1
    dy = (1 if dir else -1) * (1 if m>=0 else -1)
    np_1 = [xi + dx, yi +dy]
    np_2, np_3 = [xi, yi+dy], [xi +dx, yi] 
    if math.fabs(m) < 1:
        np_2, np_3 = np_3, np_2
    return np_1, np_2, np_3

def coef(m):
    if m >= 0 :
        result = 1 if m<1 else -1
    else:
        result = -1 if m > -1 else 1
    return result 

def _point_search(
    pi:Tuple[Union[int, float], Union[int, float]], 
    line_param: Tuple[Union[int, float], Union[int, float]], 
    to_next:int,
    allow_cross_p=False, include_pi=False):
    if type(to_next) != int:
        raise TypeError(f"\'to_next\' param must be integer, {to_next}")
    #-----------------------
    m, d = line_param
    f = f_generator(m, d)
    m_cor = coef(m)
    dir = True if to_next > 0 else False
    sgn = 1 if to_next > 0 else -1
    p = pi
    points_main =[pi] if include_pi else []
    points_sub = []
    #------------------------
    #print("to_next:", to_next)
    for i in range(0, abs(to_next)):
        p1, p2, p3 = next_points(p, m, dir)
        #print(p1, p2, p3)
        pd = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2]
        D = sgn*(m_cor) * f(*pd)
        #print("D:", D)
        if D == 1:
            points_main.append(p1)
            if allow_cross_p:
                points_sub.append(p2)
                points_sub.append(p3)
            pf = p1
        elif D == -1:
            points_main.append(p2)
            pf = p2
        else: # D == 0
            p_ = [p1[0] + (p2[0]- p[0]), p1[1] + (p2[1]- p[1])]
            points_main.append(p_)
            points_sub.append(p1)
            points_sub.append(p2)
            pf = p_ # 2 step
        p = pf
    
    #print(points_main)
    #print(points_sub)
    main = None if len(points_main) == 0 else np.array(points_main)
    sub = None if len(points_sub) == 0 else np.array(points_sub)
    return main, sub

def points(
        initial_point, 
        line_param:Tuple[float, float], 
        p_range:Tuple[int, int],
        allow_cross_p = False
        ):
        # Positive direction
        to_next = int(p_range[1])
        if to_next == 0:
            pos_main = pos_sub = None
        else:
            to_next = to_next if to_next >0 else -to_next
            pos_main, pos_sub = _point_search(
                initial_point, line_param, to_next, 
                allow_cross_p, include_pi=True) 
        # Negative direction
        to_next = p_range[0]
        if to_next == 0:
            neg_main = neg_sub = None
        else:
            to_next = to_next if to_next <0 else -to_next 
            neg_main, neg_sub = _point_search(
                initial_point, line_param, to_next, 
                allow_cross_p, include_pi=False) 

        if (pos_main is None) != (neg_main is None):
            if pos_main is None:
                main = neg_main
            else:
                main = pos_main
        elif (pos_main is not None):
            main = np.append(pos_main, neg_main, axis=0)
        else:
            main = None

        if (pos_sub is None) != (neg_sub is None):
            if pos_sub is None:
                sub = neg_sub
            else:
                sub = pos_sub
        elif (pos_sub is not None):
            sub = np.append(pos_sub, neg_sub, axis=0)
        else:
            sub = None
       
        return main, sub
