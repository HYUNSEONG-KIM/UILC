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

#
#
#

import numpy as np
import math
from scipy import optimize as op
from .hypergeo import hyper2F1

EPS = np.finfo(float).eps *1000

class Radiation:
    def lambertian(s, h, x, t):
        return h**s/(h**2 + (x-t)**2)**(s/2+1)
    def gaussian(s, h, x, t):
        return (1/(h**2 + (t-x)**2))*np.exp(- s*((t-x)/h)**2)
    def gauss_tan(s, h, x, t):
        return (1/(h**2 + (t-x)**2))*np.exp(- s*(np.tan(np.abs(t-x)/h)**2))

class Array(np.ndarray):
    def csym_index(N):
        return np.array([(i-(N-1)/2)] for i in range(0, N))
    #-------------------------------------------------
    @classmethod
    def from_arrays(cls, xarr, yarr=np.array([0.])):
        xarr = np.array(xarr)
        yarr = np.array(yarr)[::-1]

        if len(xarr.shape) != 1 or len(yarr.shape) != 1:
            raise ValueError("The given arrays, xarr and yarr, must be 1 dim array.")
        return cls(np.flip(np.array(np.meshgrid(xarr,yarr)).transpose(), axis=None))
    @classmethod
    def uniform(cls, d_t, N_t):
        d_t = tuple(d_t)
        N_t = list(N_t)

        if len(d_t) == 1:
            dx = dy =d_t[0]
        elif len(d_t) == 2:
            dx, dy = d_t
        
        if len(N_t) == 1:
            N = M = N_t[0]
        else:
            N, M = N
            N = int(N)
            M = int(M)
        
        if N == 0:
            return None
        if math.isclose(dx, 0., abs_tol =  EPS):
            return None
        
        xarr = dx*np.array([(i-(N-1)/2)] for i in range(0, N))
        yarr = dy*np.array([(j-(M-1)/2)] for j in range(0, M))
        
        #indexing = cls([[[(i-(N-1)/2), j-(M-1)/2] for i in range(0,N)] for j in range(0,M)])

        return cls.from_arrays(xarr, yarr) 
    #---------------------------------------------------
    def get_axis_list(self, axis="x"):
        shape = self.shape
        if len(shape) != 3 or len(shape) == 3 and len(self[0][0]) !=2 :
            raise ValueError("Not a 2 dimensional array of tuple.")
        N = shape[1]
        M = shape[0]
        if axis == "x" or axis == 0:
            return self[0,0:N][0:N,0]
        elif axis == "y" or axis == 1:
            return self[0:M,0][0:M,1]
        else:
            raise ValueError("axis argument must be 'x', 0 or 'y, 1 current={}".format(axis))
    def intensity_on(self, plane_points):
        pass
        return ""

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
    def get_2d_array(xarr, yarr=np.array([0.]), dim= 1):
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
    def transformMatrix(n):
        Fd = np.array([[1 if i == j else (-1 if i-j == 1 else 0) for j in range(0,n)] for i in range(0,n)])
        inFd =  np.array([[1 if i >= j else (0) for j in range(0,n)] for i in range(0,n)])

        return Fd, inFd
    def loc_to_diff(x, n):
        fd, infd = Utils.transformMatrix(n)
        d= fd.dot(x)
        w= -2*d[0]
        d0 = d[1:]
        m = math.floor(n/2) if n%2 == 0 else math.floor((n-1)/2)
        #print(m)
        return d0[0:m], w, m
    def diff_to_loc(diff, n, w):
        fd , infd = Utils.transformMatrix(n)

        if n%2 ==0:
            m = n/2
            d = np.append(diff, np.flip(diff[0:m]))
        else:
            m = (n-1)/2
            d= np.append(diff, np.flip(diff))
        np.insert(d, 0, -w/2)
        return infd.dot(d)
    def f_residual(d, **kwargs):
        W = kwargs["W"]
        h = kwargs["h"]
        n = kwargs["n"]
        xdata = kwargs["xdata"]
        ydata = kwargs["ydata"]
        m = d.size
        fd, infd = Utils.transformMatrix(n)

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
    def coefficient(s, N, M=1, shape = "L" ,approx=False):
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
    def get_nmax(s, W, H):

        xlim = W/2
        n = 2
        d = ESC.coefficient(s, n)*H
        nxe = (n-1)/2 *d

        while((nxe < xlim)): # find critical 'n' fill the given region.
            n += 1
            d = ESC.coefficient(s, n)*H
            nxe = (n-1)/2 *d
            nxm = d/2 if n%2 ==0 else d
        
        n_o = n
        # For odd and even n, their order by requiring area can be reversed.
        n1 = n-1
        d1 =  ESC.coefficient(s, n1)*H
        n1_area = (n1-1)/2 *d1
        
        n2 = n-2
        d2 = ESC.coefficient(s, n2)*H
        n2_area = (n2-1)/2 *d2
        
        n1_residual = xlim - n1_area
        n2_residual = xlim - n2_area

        n_residual = n1_residual/2 if n1_residual < n2_residual else n2_residual/2
        n = n1 if n1_residual < n2_residual else n2

        n_o_residual = math.fabs(nxe - xlim)/2 
        if n_o_residual < d:
            n = n if n_residual < n_o_residual else n_o
        return n
    
    def array(s, N, M=1, shpae="L", approx=False):
        d = ESC.coefficient(s, N, M, shpae, approx)
        xarr = d* Array.csym_index(N)
        yarr = d* Array.csym_index(M)
        return Array.from_arrays(xarr, yarr)

class OP: #Distribution optimization including bc expansion method
    def D(d, alpha, s): #Done
        return (1+ (0.5*alpha+d)**2)**(-s/2 -1) + (1+ (0.5*alpha-d)**2)**(-s/2-1) - 2*((1+d**2)**(-s/2-1))
    def Di(x, s, I0, W, H): # same with (I0/H^2 )*D(x/H, W/H, s)
        return I0/(H)**2 *OP.D(x/H, W/H, s)
    #--------------------------------------------------------
    def de(alpha, s, approx=False): #Done
        app = 0.25*alpha +  math.sqrt(2**(2/(s+2)) -1)/6
        if approx:
            return app
        else:
            r = op.root_scalar(lambda d: OP.D(d, alpha, s), x0 = app, bracket =[0, alpha/2], method = "brentq")
            return r.root
    def dm(alpha, s, de, approx=False): #Done
        app = math.sqrt(2**(2/(s+2))-1)
        if approx:
            return app
        else:
            r = op.root_scalar( lambda d: OP.D(d, alpha, s) + OP.D(alpha/2, alpha, s), x0 = app, bracket=[0,de], method="brentq")
            return r.root
    def xe(s, W, H, approx=False): #Done
        alpha = W/H
        return H*OP.de(alpha, s, approx)
    def xm(s, W, H, xe, approx=False): #Done
        if xe is None:
            xe = OP.xe(s, W, H, approx)
        de = xe/H
        alpha = W/H
        return H * OP.dm(alpha, s, de, approx)
    
    #--------------------------------------------------------
    def _bc_n_esc_max(s, H, xm, xe, status = 0): #Done, status: 0 fill Q region including P, 1: only Q
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
    #--------------------------------------------------------
    @classmethod
    def get_bc_expansion(
        pq_arr, 
        s, 
        H, 
        Wx,
        Wy=0.0,
        axis=0): # axis: (0=x), (1=y), (2= x, y all)

        shape = pq_arr.shape
        N = shape[1]
        M = shape[0]

        xarr = np.unique(Utils.get_axis_list(pq_arr, axis="x"))
        yarr = np.unique(Utils.get_axis_list(pq_arr, axis="y"))

        dx = xarr[1] - xarr[0]
        dy = yarr[1] - yarr[0] if M > 1 else 0
        
        xe = OP.xe(s, Wx, H)
        xm = OP.xm(s, Wx, H, xe)

        if (xe - xarr.max()) < 0.4*dx:
            xarr[0] = xarr[N-1] = xe\
              
        if dy != 0:
            ye = OP.xe(s, Wy, H)
            ym = OP.ym(s, Wy, H, ye)
            if (ye - yarr.max()) < dy/2:
                yarr[0] = yarr[N-1] = ye
          
        def R_region(x, s, H, W, xe, xm):
            sgn= math.copysign(1,x)
            x = math.fabs(x)
            if x> W/2 or x<0:
                raise ValueError("Argument 'x' must be in xm <= x < xe < x <= W/2 ")
            if math.isclose(x, xe, rel_tol = EPS):
                #return xe
                return None
            if math.isclose(x, W/2, rel_tol=EPS) and x<W/2:
                return xm
            if x < xm: # P region
                return W/2
            if xm < x and x < W/2: # Q region
                L = OP.D(x/H, W/H, s)
                try:
                    sol = op.root_scalar(lambda xc: OP.D(xc, W/H, s) + L, bracket=[0,W/(2*H)], method="brentq")
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
                        x_new = R_region(x,s,H,Wx,xe,xm)
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
                    y_new = R_region(y,s,H,Wy,ye,ym)
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

        return Array.from_arrays(x_arr_ex, y_arr_ex)

    def solve_linear(s, W, H, n_pre=False, k =0):
        def solve_discretized(n, s, W, H, getfd=False):
            d= W/n
            F = np.fromfunction(lambda i, j: Radiation.lambertian(s, H, d*i, d*j), (n,n), dtype=float)
            delta = np.linalg.solve(F, np.ones(n))
            if getfd:
                return delta, F
            return delta.min()

        # Get boundary value of positiive solution
        napp = math.floor(3/ hyper2F1(1/2, (s+2)/2, 3/2, - (W/(2*H))**2))
        state = 0
        termination = False
        while(not termination):
            minvalue = solve_discretized(n, s, W, H)
            if minvalue <0:
                if state ==2:
                    napp -= 1
                    termination = True
                else:
                    napp -=1
                    state =1
            else:
                if state ==1:
                    termination =True
                else:
                    napp +=1
                    state=2
        
        dimension_n = n_pre if n_pre else napp - k
        
        if dimension_n >= napp:
            # Solve discretized
            n = dimension_n
            d= W/n
            delta, F = solve_discretized(n, s, W, H, getfd=True)
            position =np.array([-W/2+d/2+(i-1)*d for i in range(1,n+1)])
            return delta, position, F
        else:
            # using active set method and solve NNLS

        

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
        d, w, m = Utils.loc_to_diff(Utils.get_axis_list(arr), n)
        fd, infd = Utils.transformMatrix(n)
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

