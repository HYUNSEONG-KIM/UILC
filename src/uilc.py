#Copyright (C) 2022 Hyun Seong Kim (김현성)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
import math
from scipy import special as sp
from scipy import optimize as op
from hypergeo import hyper2F1

EPS = np.finfo(float).eps *1000


class Utils: # basic Utils and functions including mathematical routines
    def uniformarray(dx, N, dy=0.0, M=1, diff = False): #return uniformly distributed array of which distance between elements is 'd'
        indexing = np.array([[[(i-(N-1)/2), j-(M-1)/2] for i in range(0,N)] for j in range(0,M)])
        if diff:
            for i in range(0,N):
                indexing[0:M,i][0:M,0] *= dx
                indexing[0:M,i][0:M,1] *= dy
            return indexing
        else:
            return np.multiply(indexing, dx)
    def get_axis_list(arr, axis="x"): 
        shape = arr.shape
        N = shape[1]
        M = shape[0]
        if axis == "x":
            return arr[0,0:N][0:N,0]
        elif axis == "y":
            return arr[0:M,0][0:M,1]
        else:
            raise ValueError("axis argument must be 'x' or 'y, current={}".format(axis))
    def get_2d_array(xarr, yarr, dim= 1):
        if dim ==2:
            xarr = Utils.get_axis_list(xarr)
            yarr = Utils.get_axis_list(yarr)

        return  np.flip(np.array(np.meshgrid(xarr,yarr)).transpose(), axis=None)
    def getlocation(i, j, array): # (x-order, y-order) conversion for matrix row, column order indexing.
        return array[j,i]
    def gauss_distribution(x, y, arr, h): #approximating Dirac Delta function
        M = arr.shape[0]
        N = arr.shape[1]
        result = 0.0
        for i in range(0,N):
            for j in range(0,M):
                d = (x-arr[j,i][0])**2 + (y - arr[j,i][1])**2
                result += (1/(math.fabs(h)*np.sqrt(2*math.pi)))*np.exp(-d/(2*h**2))
        return result
    def dirac(x, y, arr):
        M = arr.shape[0]
        N = arr.shape[1]
        result = 0.0
        for i in range(0,N):
            for j in range(0,M):
                d = 1 if (x-arr[j,i][0])**2 < EPS and (y - arr[j,i][1])**2 < EPS else 0
                result += d
        return result
    def intensity_function(s, h, x, t):
        return h**s/(h**2 + (x-t)**2)**(s/2+1)
    def lambertian_distribution(x, y, arr, s, h):
        M = arr.shape[0]
        N = arr.shape[1]
        result =0.0
        for i in range(0,N):
            for j in range(0,M):
                d = (x-arr[j,i][0])**2 + (y - arr[j,i][1])**2

                result += h**s / (h**2 + d)**(s/2 +1)
        return result
    def evaluation_factors(arr):
        stdE = (arr.max() - arr.min())/arr.mean() # |E_max -E_min|/E_mean
        rmse = (arr.std())/arr.mean() # std/mean
        return stdE, rmse
#=============================================================================================================================================
    def ls_transformMatrix(n):
        Fd = np.array([[1 if i == j else (-1 if i-j == 1 else 0) for j in range(0,n)] for i in range(0,n)])
        inFd =  np.array([[1 if i >= j else (0) for j in range(0,n)] for i in range(0,n)])

        return Fd, inFd
    def ls_loc_to_diff(x, n):
        fd, infd = Utils.ls_transformMatrix(n)
        d= fd.dot(x)
        w= -2*d[0]
        d0 = d[1:]
        m = math.floor(n/2) if n%2 == 0 else math.floor((n-1)/2)
        #print(m)
        return d0[0:m], w, m
    def ls_diff_to_loc(diff, n, w):
        fd , infd = Utils.ls_transformMatrix(n)

        if n%2 ==0:
            m = n/2
            d = np.append(diff, np.flip(diff[0:m]))
        else:
            m = (n-1)/2
            d= np.append(diff, np.flip(diff))
        np.insert(d, 0, -w/2)
        return infd.dot(d)
    def ls_f_residual(d, **kwargs):
        W = kwargs["W"]
        h = kwargs["h"]
        n = kwargs["n"]
        xdata = kwargs["xdata"]
        ydata = kwargs["ydata"]
        m = d.size
        fd, infd = Utils.ls_transformMatrix(n)

        if n%2 ==0:
            d0 = np.append(d[0:m], np.flip(d[0:m-1]))
            pass
        else:
            d0 = np.append(d,np.flip(d))

        d0 = np.insert(d0, 0, -W/2)
        xi = infd.dot(d0)

        location = np.array([[[x, 0] for x in xi]])

        return (ydata - Utils.gauss_distribution(xdata, np.array([[[0]]]), location,h))[0][0]


class ESC: #Expanded Sparrow Criterion
    def __init__(self):
        pass
    @classmethod
    def coefficient(cls, s, N, M=1, shape = "L" ,approx=False):

        #Value check : s>= 1.0, N, M are natural numbers, shape ="L" or "R"
        if not isinstance(s,(float, int)) or not isinstance(N, int) or not isinstance(M, int):
            message = 'Check the types of arguments s = int or float >=1, N and M = int>0 \n Current types are s={}, N ={}, M={}'.format(type(s), type(N), type(M))
            raise TypeError(message)
        if s <1.0 or N <0 or M <0:
            message = 'Domain Error: s >= 1, N and M >=1 \n Current arguments: s={}, N ={}, M={}'.format((s), (N), (M))
            raise ValueError(message)


        def linearfunction(D, N):
            y =0.0
            for i in range(1,N+1):
                y += (1-(s+3)*(N+1-2*i)**2 * (D**2)/4)*((N+1-2*i)**2 * (D**2)/4 +1)**(-(s+6)/2)
            return y
        def rectangularfunction(D, N, M):
            y =0.0
            for i in range(1,N+1):
                for j in range(1, M+1):
                    y += (((N+1-2*i)**2 + (M+1-2*j)**2)*(D**2/4)+1)**(-(s/2+3.0)) * (1-((s+3)*(N+1-2*i)**2 -(M+1-2*j)**2)*(D**2)/4)
            return y
        
        if shape == "L":
            if N==2:
                return math.sqrt(4/(s+3))
            if approx ==True and (N>4 and s>30):
                return math.sqrt(3.2773/(s+4.2539))
            
            if N%2 ==0: #Even number of LED
                sol = op.root_scalar(lambda D: linearfunction(D ,N), bracket=[0,1], method = "brentq")
                return sol.root
            else: #Odd number of LED
                res = op.minimize_scalar(lambda D: linearfunction(D, N), bounds=(0,1), method = "bounded")
                return  res.x


        if shape == "R":
            if N==2 and N == M:
                return math.sqrt(4/(s+2))
            if approx == True and (N > 4 and M > 4 and s>30):
                return math.sqrt(1.2125/(s-3.349))
            try:
                sol = op.root_scalar(lambda D: rectangularfunction(D, N, M), bracket=[0,1],method="brentq")
                if sol.converged == False:
                    res = op.minimize_scalar(lambda D: rectangularfunction(D, N, M), bounds=(0,1), method = "bounded")
                    return res.x
                else:
                    return sol.root
            except:
                res = op.minimize_scalar(lambda D: rectangularfunction(D, N, M), bounds=(0,1), method = "bounded")
                return res.x
        
        
        raise ValueError("\"shape\" argument must be \"L\" or \"R\" current value is {}".format(shape))

    @classmethod
    def get_nmax(cls, s, W, H):
        xlim = W/2
        n = 2
        d = cls.coefficient(s, n)*H
        nxe = (n-1)/2 *d

        while((nxe < xlim)):
            n += 1
            d = cls.coefficient(s, n)*H
            nxe = (n-1)/2 *d
            nxm = d/2 if n%2 ==0 else d

        n1 = n-1
        d1 =  cls.coefficient(s, n1)*H
        nxe1 = (n1-1)/2 *d1
        
        n2 = n-2
        d2 = cls.coefficient(s, n2)*H
        nxe2 = (n2-1)/2 *d2
        
        di1 = xlim - nxe1
        di2 = xlim - nxe2

        n = n1 if di1 < di2 else n2

        return n
class DISOP: #Distribution optimization including bc expansion method
    def __init__(self, s, W, H):
        pass
    @classmethod
    def D(cls, d, alpha, s): #Done
        return (1+ (0.5*alpha+d)**2)**(-s/2 -1) + (1+ (0.5*alpha-d)**2)**(-s/2-1) - 2*((1+d**2)**(-s/2-1))
    @classmethod
    def find_de(cls, alpha, s, approx=False): #Done
        app = 0.25*alpha +  math.sqrt(2**(2/(s+2)) -1)/6
        if approx:
            return app
        else:
            r = op.root_scalar(lambda d: cls.D(d, alpha, s), x0 = app, bracket =[0, alpha/2], method = "brentq")
            return r.root
    @classmethod
    def find_dm(cls, alpha, s, de, approx=False): #Done
        app = math.sqrt(2**(2/(s+2))-1)
        if approx:
            return app
        else:
            r = op.root_scalar( lambda d: cls.D(d, alpha, s) + cls.D(alpha/2, alpha, s), x0 = app, bracket=[0,de], method="brentq")
            return r.root
    @classmethod
    def find_xe(cls, W, H, s, approx=False): #Done
        alpha = W/H
        return cls.find_de(alpha, s, approx) *H
    @classmethod
    def find_xm(cls, W, H, s, xe, approx=False): #Done
        de = xe/H
        alpha = W/H
        return cls.find_dm(alpha, s, de, approx) *H
    @classmethod
    def Di(cls, s, I0, x, W, H): # same with (I0/H^2 )*D(x/H, W/H, s)
        return I0/(H)**2 *cls.D(x/H, W/H, s)
    @classmethod
    def get_n_ESC_max(cls, s, H, xm, xe, status = 0): #Done, status: 0 fill Q region including P, 1: only Q
        n =2
        d = ESC.coefficient(s, n)*H
        nxe = (n-1)/2 *d
        nxm = d/2 if n%2 ==0 else d
        while((nxe < xe-d/2)):
            n += 1
            d = ESC.coefficient(s, n)*H
            nxe = (n-1)/2 *d
            nxm = d/2 if n%2 ==0 else d
            if status ==1 and nxm > xm:
                break
            
        if nxe > xe and n > 2:
                n1 = n-1
                d1 =  ESC.coefficient(s, n1)*H
                nxe1 = (n1-1)/2 *d1
                n2 = n-2
                d2 = ESC.coefficient(s, n2)*H
                nxe2 = (n2-1)/2 *d2
                di1 = xe - nxe1
                di2 = xe - nxe2
                n = n1 if di1 < di2 else n2

        return n

#=============================================================================================
    @classmethod
    def get_bc_expansion(cls, ESC_arr, s, H, Wx, xe, xm, Wy=0.0, ye=0.0, ym=0.0, axis=0): # axis: (0=x), (1=y), (2= x, y all)
        shape = ESC_arr.shape
        N = shape[1]
        M = shape[0]

        xarr = np.unique(Utils.get_axis_list(ESC_arr, axis="x"))
        yarr = np.unique(Utils.get_axis_list(ESC_arr, axis="y"))

        dx = xarr[1] - xarr[0]
        dy = yarr[1] - yarr[0] if M > 1 else 0

        if (xe - xarr.max()) < 2*dx/5:
            xarr[0] = xarr[N-1] = xe
        if dy != 0 and (ye - yarr.max()) < dy/2:
            yarr[0] = yarr[N-1] = ye

        def core(x, s, H, W, xe, xm):
            sgn= math.copysign(1,x)
            x = math.fabs(x)
            if x> W/2 or x<0:
                raise ValueError("Argument 'x' must be in xm <= x < xe < x <= W/2 ")
            if math.isclose(x, xe, rel_tol = EPS):
                return xe
            if math.isclose(x, W/2, rel_tol=EPS) and x<W/2:
                return xm
            if x < xm:
                return W/2
            if xm < x and x < W/2:
                L = cls.D(x/H, W/H, s)
                try:
                    sol = op.root_scalar(lambda xc: cls.D(xc, W/H, s) + L, bracket=[0,W/(2*H)], method="brentq")
                except ValueError:
                    raise ValueError("x: {}, xe:{}, xm:{}".format(x, xe, xm))
                return sol.root * H * sgn
        x_ex = []
        y_ex = []
        if axis == 0 or axis ==2:
            for i, x in enumerate(xarr):
                    if x > xe:
                        np.delete(xarr,i)
                    else:
                        x_new = core(x,s,H,Wx,xe,xm)
                        if math.fabs(x_new - x) < dx/2:
                            np.delete(xarr,i)
                            x_ex.append(xe)
                        else:
                            x_ex.append(x_new)
                    
        if axis == 1 or axis == 2:
            for i, y in enumerate(yarr):
                if y > ye:
                    np.delete(yarr,i)
                else:
                    y_new = core(y,s,H,Wy,ye,ym)
                    if math.fabs(y_new - y) < dy/3:
                        np.delete(yarr, i)
                        y_ex.append(ye)
                    else:
                        y_ex.append(y_new)
        
        if len(y_ex) == 0:
            y_ex= [0]

        x_ex =np.unique(np.array(x_ex))
        x_ex.sort()
        y_ex = np.unique(np.array(y_ex))
        y_ex.sort()

        x_arr_ex = np.append(xarr, x_ex)
        y_arr_ex = np.append(yarr, y_ex)

        x_arr_ex = np.unique(np.append(-x_arr_ex,x_arr_ex))
        y_arr_ex = np.unique(np.append(-y_arr_ex,y_arr_ex))

        x_arr_ex.sort()
        y_arr_ex.sort()

        return np.array([[[x,y] for x in x_arr_ex] for y in y_arr_ex])

    @classmethod
    def solve_linear(cls, s, W, H, n_pre=False, k =0):
        def solve_discretized(n, s, W, H, getfd=False):
            d= W/n
            F = np.fromfunction(lambda i, j: Utils.intensity_function(s, H, d*i, d*j), (n,n), dtype=float)
            delta = np.linalg.solve(F, np.ones(n))
            if getfd:
                return delta, F
            return delta.min()

        if n_pre:
            n = n_pre
        else:
            napp = 3/ Utils.hyper2F1(1/2, (s+2)/2, 3/2, - (W/(2*H))**2)
            n = math.floor(napp)
            state = 0
            termination = False
            while(not termination):
                minvalue = solve_discretized(n, s, W, H)
                if minvalue <0:
                    if state ==2:
                        n -= 1
                        termination = True
                    else:
                        n-=1
                        state =1
                else:
                    if state ==1:
                        termination =True
                    else:
                        n+=1
                        state=2
            n = n-k
        
        d= W/n
        delta, F = solve_discretized(n, s, W, H, getfd=True)
        position =np.array([-W/2+d/2+(i-1)*d for i in range(1,n+1)])
        return delta, position, F

    @classmethod
    def solve_nnls(cls, s, W, H, n_nnls = 1, mean=True): #linear, nnls
        d= W/n_nnls
        F = np.fromfunction(lambda i, j: Utils.intensity_function(s, H, d*i, d*j), (n_nnls,n_nnls), dtype=float)
        delta = op.nnls(F,np.ones(n_nnls))[0]
        position =np.array([-W/2+d/2+(i-1)*d for i in range(1,n_nnls+1)])
        #Get meaningful points
        if mean:
            therhold = 0.01 * delta.max()
            #therhold = 2* delta.mean()
            position = position[np.argwhere(delta>therhold )[0:,0]]
            delta = delta[delta > therhold ]
        return delta, position, F
    @classmethod
    def nomarlization_lq(cls, arr, xdata, ydata, n, h=False, W=False):
        d, w, m = Utils.ls_loc_to_diff(Utils.get_axis_list(arr), n)
        fd, infd = Utils.ls_transformMatrix(n)
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

        sol = op.least_squares(Utils.ls_f_residual, x0 =bc_x0, kwargs=kwargs)
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

