
# Dot searching algrithm using Bresenham line algorithm
# refer:
#   Bresenham line algorithm
#       - https://www.cs.helsinki.fi/group/goa/mallinnus/lines/bresenh.html
#       - https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm, 10, Feb, 2023

import math
from typing import Tuple, Union

import numpy as np

from uilc.utils import float_eps

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
        raise TypeError("\'to_next\' param must be integer")
    #-----------------------
    m, d = line_param
    f = f_generator(m, d)
    m_cor = coef(m)
    dir = True if to_next > 0 else False
    p = pi
    points_main =[pi] if include_pi else []
    points_sub = []
    #------------------------
    for i in range(0, to_next):
        p1, p2, p3 = next_points(p, m, dir)
        pd = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2]
        D = (m_cor) * f(*pd)
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
    return np.array(points_main), np.array(points_sub)

def points(
        initial_point, 
        line_param:Tuple[float, float], 
        p_range:Tuple[int, int],
        allow_cross_p = False
        ):
        to_next = p_range[0]
        if to_next == 0:
            pos_main = pos_sub = np.array([])
        else:
            to_next = to_next if to_next >0 else -to_next
            pos_main, pos_sub = _point_search(
                initial_point, line_param, to_next, 
                allow_cross_p, include_pi=True) 
        to_next = p_range[0]
        if to_next == 0:
            neg_main = neg_sub = np.array([])
        else:
            to_next = to_next if to_next <0 else -to_next 
            neg_main, neg_sub = _point_search(
                initial_point, line_param, to_next, 
                allow_cross_p, include_pi=False) 
        return pos_main + neg_main, pos_sub + neg_sub

#-----------------------------------------------#
# Visualization                                 #
#-----------------------------------------------#

# from matplotlib import pyplot as plt
# pi = [0, 0] # center point
# ms = [ 0, 0.5, 1, 2, -0.6, -1, -3]
# bs = [pi[1] - m*pi[0] for m in ms ]

# n =30
# x = 1.2*n
# line_points = [[[0, x], [line(0, m, b), line(x, m, b)]] for m, b in zip(ms, bs)]
# p_range = (30, 30)
# int_points = [points(pi, line_param=[m, b], p_range=p_range) for m, b in zip(ms, bs)]
# line_pixel_points = [[p[0][:-1, 0], p[0][:-1, 1]] for p in int_points]
#
# fig = plt.figure(figsize=(8,8))
# ax = fig.add_subplot()
# for l_p in line_points:
#     ax.plot(l_p[0], l_p[1])
# for point in line_pixel_point:
#     ax.scatter(point[0], point[1])
# ax.set_xlim(-(n/2+0.5), n/2+0.5)
# ax.set_ylim(-(n/2+0.5), n/2+0.5)
# ax.grid()
# plt.show()