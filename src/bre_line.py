import math
from typing import Tuple
from copy import deepcopy

import numpy as np


EPS = np.finfo(float).eps

def line(x, m, b):
    return m * x + b

def plane_deter(point, m, d, dx = 1):
    a = m* dx
    b = - dx
    c = dx* d
    x, y = point
    result = (a*x + b*y + c)
    if math.fabs(result) < EPS:
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
        if m < 1:
            result = 1
        else:
            result = -1
    else:
        if m > -1:
            result = -1
        else:
            result = 1
    return result 


def points(
        initial_point, 
        line_param:Tuple[float, float], 
        p_range:Tuple[int, int],
        allow_cross_thick = False
        ):
        m, d = line_param
        f = f_generator(m, d)
        m_cor = coef(m)
        points_main =[initial_point]
        points_sub = []
        
        # positive direction
        dir = True
        p = deepcopy(initial_point)
        for i in range(0, p_range[1]):
            p1, p2, p3 = next_points(p, m, dir)
            pd = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2]
            D = (m_cor) * f(*pd)
            if D == 1:
                points_main.append(p1)
                if allow_cross_thick:
                    points_main.append(p2)
                    points_main.append(p3)
                pf = p1
            elif D == -1:
                points_main.append(p2)
                pf = p2
            else:
                points_main.append(p1)
                points_main.append(p2)
                p_ = [p1[0] + (p2[0]- p[0]), p1[1] + (p2[1]- p[1])]
                points_main.append(p_)
                pf = p_ # 2 step
            p = pf
        # negative direction
        dir = False
        p = deepcopy(initial_point)
        for i in range(0, p_range[0]):
            p1, p2, p3 = next_points(p, m, dir)
            pd = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2]
            D = (m_cor) * ( -f(*pd)) 
            if D == 1:
                points_main.append(p1)
                if allow_cross_thick:
                    points_main.append(p2)
                    points_main.append(p3)
                pf = p1
            elif D == -1:
                points_main.append(p2)
                pf = p2
            else:
                points_main.append(p1)
                points_main.append(p2)
                p_ = [p1[0] + (p2[0]- p[0]), p1[1] + (p2[1]- p[1])]
                points_main.append(p_)
                pf = p_ # 2 step
                
            p = pf
        return np.array(points_main)

#-------------------------------------------------------------------
# 
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