
import math

import numpy as np
from scipy import optimize as op

from uilc import PositionArray
from uilc.utils import float_eps

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

def _r_region_x(x, s, H, W, xe, xm):
    sgn= math.copysign(1,x)
    x = math.fabs(x)
    if x> W/2 or x<0:
        raise ValueError("Argument 'x' must be in xm <= x < xe < x <= W/2 ")
    if math.isclose(x, xe, abs_tol = float_eps):
        return None
    if math.isclose(x, W/2, abs_tol=float_eps) and x<W/2:
        return xm
    if x < xm: # P region
        return W/2
    if xm < x and x < W/2: # Q region
        L = _D_f(x/H, W/H, s)
        try:
            sol = op.root_scalar(lambda xc: _D_f(xc, W/H, s) + L, bracket=[0,W/(2*H)], method="brentq")
        except ValueError:
            raise ValueError("x: {}, xe:{}, xm:{}".format(x, xe, xm))
        return sol.root * H * sgn

def _get_r_points(xarr, dx, s, H, W, xe, xm, thershold=0.7):
    x_ex = []
    for i, x in enumerate(xarr):
        if x > xe:
            np.delete(xarr,i)
        else:
            x_new = OP._r_region_x(x,s,H,W,xe,xm)
            if x_new is None:
                xarr[i] = xe
            else:
                if math.fabs(x_new - x) < dx*thershold:
                    np.delete(xarr,i)
                    x_ex.append(xe)
                else:
                    x_ex.append(x_new)
    return x_ex


def fill_rq(s, H, xm, xe, status = 0): 
    #Done, status: 0 fill Q region including P, 1: only Q
    n =2
    d = ESC.coefficient(s, n)[0]*H
    nxe = (n-1)/2 *d
    nxm = d/2 if n%2 ==0 else d
    while((nxe < xe-d/2)):
        n += 1
        d = ESC.coefficient(s, n)[0]*H
        nxe = (n-1)/2 *d
        nxm = d/2 if n%2 ==0 else d
        if status ==1 and nxm > xm:
            break
        
    if nxe > xe and n > 2:
            n1 = n-1
            d1 =  ESC.coefficient(s, n1)[0]*H
            nxe1 = (n1-1)/2 *d1
            n2 = n-2
            d2 = ESC.coefficient(s, n2)[0]*H
            nxe2 = (n2-1)/2 *d2
            di1 = xe - nxe1
            di2 = xe - nxe2
            n = n1 if di1 < di2 else n2
    return n


def get_bc_expansion(
        pq_arr:PositionArray,
        s:float, 
        H:float, 
        Wx:float,
        Wy:float,
        expand=False):

        shape = pq_arr.shape
        N = shape[0]
        M = shape[1]

        xarr = np.unique(PositionArray.get_axis_list(pq_arr, axis="x"))
        yarr = np.unique(PositionArray.get_axis_list(pq_arr, axis="y"))

        dx = xarr[1] - xarr[0]
        dy = yarr[1] - yarr[0]
        
        xe = get_xe(s, Wx, H)
        xm = get_xm(s, Wx, H, xe)
        ye = get_xe(s, Wx, H)
        ym = get_xm(s, Wx, H, ye)

        if (xe - xarr.max()) >= 0.4*dx:
            xarr[0] = xarr[N-1] = xe
        if (ye - yarr.max()) >= 0.4*dy:
            yarr[0] = yarr[M-1] = ye
        
        x_ex = _get_r_points(xarr, dx, s, H, Wx, xe, xm)
        y_ex = _get_r_points(yarr, dy, s, H, Wy, ye, ym)

        #print(f"x_ex:{x_ex}")
        #print(f"xarr:{xarr}")
        #print(f"dx:{dx}")
        #print(f"y_ex:{y_ex}")
        #print(f"yarr:{yarr}")
        #print(f"dy:{dy}")
        #-----Do below code
        x_ex =np.unique(np.array(x_ex))
        x_ex.sort()
        x_arr_ex = np.append(xarr, x_ex)
        x_arr_ex = np.unique(np.append(-x_arr_ex,x_arr_ex))
        x_arr_ex.sort()

        y_ex =np.unique(np.array(y_ex))
        y_ex.sort()
        y_arr_ex = np.append(yarr, y_ex)
        y_arr_ex = np.unique(np.append(-y_arr_ex,y_arr_ex))
        y_arr_ex.sort()

        return x_arr_ex, y_arr_ex