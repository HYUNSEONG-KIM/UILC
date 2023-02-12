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

import math
from copy import deepcopy
from collections.abc import Iterable
from typing import Tuple

import numpy as np
from scipy import optimize as op

from uilc.utils.hypergeo import hyper2F1
from uilc.utils.bresenham import points 
from uilc.utils import float_eps, d2, csym_index

class Radiation:
    def lambertian(s, h, d, inv=True):
        r = (h**2 + d)**-1 if inv else 1
        return r/(1 + d/(h**2))**(s/2)
    def gaussian(s, h, d, inv = False):
        r =  (1/(h**2 + d)) if inv else 1.
        return r*np.exp(- s*(np.sqrt(d)/h)**2)
    def gauss_tan(s, h, d, inv=False):
        r =  (1/(h**2 + d)) if inv else 1.
        return r*np.exp(- s*(np.tan(np.abs(np.sqrt(d))/h)**2))
    def dirac(d, ep=1000*float_eps):
        return ((np.sqrt(np.pi)*ep))**(-1) * np.exp(- d/(ep)**2)

class PositionArray(np.ndarray):
    def __new__(cls, input_array, *args, **kwargs):
        obj = np.asarray(input_array, *args, **kwargs).view(cls)
        return obj
    def __array_wrap__(self, out_arr, context=None):
        return super().__array_wrap__(self, out_arr, context)
    def __array_ufunc__(self, ufunc, method, *inputs, out=None, **kwargs):
        args = []
        in_no = []
        for i, input_ in enumerate(inputs):
            if isinstance(input_, PositionArray):
                in_no.append(i)
                args.append(input_.view(np.ndarray))
            else:
                args.append(input_)

        outputs = out
        out_no = []
        if outputs:
            out_args = []
            for j, output in enumerate(outputs):
                if isinstance(output, PositionArray):
                    out_no.append(j)
                    out_args.append(output.view(np.ndarray))
                else:
                    out_args.append(output)
            kwargs['out'] = tuple(out_args)
        else:
            outputs = (None,) * ufunc.nout

        info = {}
        if in_no:
            info['inputs'] = in_no
        if out_no:
            info['outputs'] = out_no

        results = super().__array_ufunc__(ufunc, method, *args, **kwargs)
        if results is NotImplemented:
            return NotImplemented

        if method == 'at':
            if isinstance(inputs[0], PositionArray):
                inputs[0].info = info
            return

        if ufunc.nout == 1:
            results = (results,)

        results = tuple((np.asarray(result).view(PositionArray)
                         if output is None else output)
                        for result, output in zip(results, outputs))
        if results and isinstance(results[0], PositionArray):
            results[0].info = info

        return results[0] if len(results) == 1 else results
    #-------------------------------------------------
    @classmethod
    def from_arrays(cls, xarr, yarr=np.array([0.])):
        xarr = np.array(xarr)[::-1]
        yarr = np.array(yarr)

        if len(xarr.shape) != 1 or len(yarr.shape) != 1:
            raise ValueError("The given arrays, xarr and yarr, must be 1 dim array.")
        return cls(np.flip(np.array(np.meshgrid(yarr,xarr)).transpose(), axis=None))
    @classmethod
    def from_meshgrid(cls, meshx, meshy=None, indexing="xy"):
        if meshy is None and len(meshx) ==2:
            meshx, meshy = meshx
        if meshx.shape != meshy.shape:
            raise ValueError("Two mesh do not share dimension.")
        if indexing == "ij":
            meshx = meshx.transpose()
            meshy = meshy.transpose()
        elif indexing == "xy":
            pass
        else:
            raise ValueError("Indexing must be \'ij\' or \'xy\'.")
        mesh = (meshy, meshx)
        arr = np.flip(np.array(mesh).transpose(), axis=None)[::-1]

        return cls(np.transpose(arr, axes=(1, 0, 2)))
    @classmethod
    def uniform(cls, d_t, N_t):
        if isinstance(d_t, Iterable):
            d_t = list(d_t)
        else:
            d_t = [d_t, d_t]
        if isinstance(N_t, Iterable):
            N_t = list(N_t)
        else:
            N_t = [N_t, N_t]

        if len(d_t) == 1:
            dx = dy =d_t[0]
        elif len(d_t) == 2:
            dx, dy = d_t
        
        if len(N_t) == 1:
            N = M = N_t[0]
        else:
            N, M = N_t
            N = int(N)
            M = int(M)
        
        if N == 0:
            return None
        if math.isclose(dx, 0., abs_tol =  float_eps):
            return None
        
        xarr = dx*(cls.csym_index(N))
        yarr = dy*(cls.csym_index(M))
        
        #indexing = cls([[[(i-(N-1)/2), j-(M-1)/2] for i in range(0,N)] for j in range(0,M)])

        return cls.from_arrays(xarr, yarr) 
    @classmethod
    def uniform_fill(cls, W, N_t):
        Nx, Ny = N_t
        Wx, Wy = W

        dx = Wx/(Nx-1)
        dy = Wy/(Ny-1)

        return cls.uniform((dx, dy), N_t)
    #---------------------------------------------------
    def check_dim(self):
        shape = self.shape
        if len(shape) != 3 or len(shape) == 3 and len(self[0][0]) !=2 :
            raise ValueError("Not a 2 dimensional array of tuple.")
    def get_axis_list(self, axis="x"):
        self.check_dim()
        
        N = self.shape[1]
        M = self.shape[0]
        if axis == "x" or axis == 0:
            return self[0,0:N][0:N,0]
        elif axis == "y" or axis == 1:
            return self[0:M,0][0:M,1]
        elif axis == "z" or axis == 2:
            return self[0,0:N][0:N,0], self[0:M,0][0:M,1]
        else:
            raise ValueError("axis argument must be 'x', 0 or 'y, 1 current={}".format(axis))
    def to_meshgrid(self, indexing = "xy"):
        xarr = self.get_axis_list(axis = "x")
        yarr = self.get_axis_list(axis = "y")
        return np.meshgrid(xarr, yarr, indexing=indexing)
    def intensity_on(self, plane_meshes:np.ndarray, radiation_pattern:callable):
        self.check_dim()
        X, Y = plane_meshes
        result = np.zeros(shape= X.shape)
        for line in self:
            for point in line:
                p_x, p_y = point
                d = (X- p_x)**2 + (Y-p_y)**2
                result += radiation_pattern(d)
        return result

    @property
    def area(self):
        xarr = self.get_axis_list(axis = "x")
        yarr = self.get_axis_list(axis = "y")
        return (xarr.max()-xarr.min(), yarr.max()-yarr.min())

class ESC: #Expanded Sparrow Criterion
    def _linear(D, s, N):
        y =0.0
        for i in range(1,N+1):
            y += (1-(s+3)*(N+1-2*i)**2 * (D**2)/4)*((N+1-2*i)**2 * (D**2)/4 +1)**(-(s+6)/2)
        return y
    def _rectangular(D, s, N, M):
        y =0.0
        for i in range(1,N+1):
            for j in range(1, M+1):
                y += (((N+1-2*i)**2 + (M+1-2*j)**2)*(D**2/4)+1)**(-(s/2+3.0)) * (1-((s+3)*(N+1-2*i)**2 -(M+1-2*j)**2)*(D**2)/4)
        return y
    def _coefficient_linear(s, N, approx=False):
        if N == 2:
            cof = math.sqrt(4/(s+3))
        elif approx and (N >4 and s >30):
            cof = math.sqrt(3.2773/(s+4.2539))
        elif N%2 == 0 :
            sol = op.root_scalar(lambda D: ESC._linear(D, s, N), bracket=[0,1], method = "brentq")
            cof = sol.root
        else:
            res = op.minimize_scalar(lambda D: ESC._linear(D, s, N), bounds=(0,1), method = "bounded")
            cof = res.x
        return cof
    def _coefficient_rectangular(s, N, M, approx=False):
        if M > N:
            N, M = M, N
        if N==2 and N == M:
            cof= math.sqrt(4/(s+2))
        if approx == True and (N > 4 and M > 4 and s>30):
            cof= math.sqrt(1.2125/(s-3.349))
        else:
            try:
                sol = op.root_scalar(lambda D: ESC._rectangular(D, s, N, M), bracket=[0,1],method="brentq")
                if sol.converged == False:
                    res = op.minimize_scalar(lambda D: ESC._rectangular(D, s, N, M), bounds=(0,1), method = "bounded")
                    cof= res.x
                else:
                    cof= sol.root
            except:
                res = op.minimize_scalar(lambda D: ESC._rectangular(D, s, N, M), bounds=(0,1), method = "bounded")
                cof= res.x
        return cof
    def coefficient(
        s, 
        N, M=1, 
        shape = "L",
        approx=False) -> Tuple[float, float]:
        #Value check : s>= 1.0, N, M are natural numbers, shape ="L" or "R"
        #if not isinstance(s,(float, int)) or not isinstance(N, (int, np.int32)) or not isinstance(M, (int, np.int32)):
        #    message = 'Check the types of arguments s = int or float >=1, N and M = int>0 \n Current types are s={}, N ={}, M={}'.format(type(s), type(N), type(M))
        #    raise TypeError(message)
        if s <1.0 or N <0 or M <0:
            message = 'Domain Error: s >= 1, N and M >=1 \n Current arguments: s={}, N ={}, M={}'.format((s), (N), (M))
            raise ValueError(message)
        
        if shape == "L":
            cof_x = ESC._coefficient_linear(s, N, approx)
            cof_y = None if M == 1 else ESC._coefficient_linear(s, M, approx)
            return cof_x, cof_y
        if shape == "R":
            cof_x = cof_y = ESC._coefficient_rectangular(s, N, M, approx)
            return cof_x, cof_y

        raise ValueError("\"shape\" argument must be \"L\" or \"R\" current value is {}".format(shape))

    def _linear_nmax(s, W, H, thershold):
        W = W[0]
        xlim = W/2
        n = 2
        d = ESC.coefficient(s, n)[0]*H
        nxe = (n-1)/2 *d

        while((nxe < xlim)): # find critical 'n' fill the given region.
            n += 1
            d = ESC.coefficient(s, n)[0]*H
            nxe = (n-1)/2 *d
            #nxm = d/2 if n%2 ==0 else d

        n_o = n
        # For odd and even n, their order by requiring area can be reversed.
        # Below codes are choosing process by the given thershold value.
        n1 = n-1
        d1 =  ESC.coefficient(s, n1)[0] *H
        n1_area = (n1-1)/2 *d1

        n2 = n-2
        d2 = ESC.coefficient(s, n2)[0]*H
        n2_area = (n2-1)/2 *d2

        n1_residual = xlim - n1_area # By the loop termination condition residuals are guaranted to be positive.
        n2_residual = xlim - n2_area
        n_o_residual = math.fabs(nxe - xlim)
        #--------------------------------------------------------------------------
        n_index = [n1, n2, n_o]
        resi = [
            n1_residual,
            n2_residual,
            n_o_residual
        ]
        thershold_conditions = [
            n1_residual < d1 * thershold,
            n2_residual < d2 * thershold,
            n_o_residual < d * thershold
        ]

        thers_index = []
        for i, thershold_value in enumerate(thershold_conditions):
            if thershold_value:
                thers_index.append(i)
        if len(thers_index) == 0:
            ni, resi = (n1, n1_residual) if n1_residual < n2_residual else (n2, n2_residual)
            n = n_o if n_o_residual < resi else ni
            thers_condition = False
        elif len(thers_index) >1:
            thers_condition = True
            j = thers_index[0]
            n_j = n_index[j]
            n_residual = resi[j]
            for i in thers_index:
                if i == j:
                    pass
                else:
                    if resi[i] < n_residual:
                        n_j = n_index[i]
                        n_residual = resi[i]
            n = n_j                
        else:
            thers_condition = True
            n = n_index[thers_index[0]]

        return (n, n, 1, thers_condition) # (nmax, n_x, n_y, thershold value satisfaction), n_x * n_y = nmax
    
    def _area(s, nx, ny=1, shape = "L"):
        dx, dy = ESC.coefficient(s, nx, ny)
        dy = 0 if dy is None else dy
        if shape == "L":
            return (nx-1)*dx, 0
        elif shape == "R" and ny >=1:
            return (nx -1)*dx, (ny-1)*dy
        else:
            raise ValueError(f"shape:{shape}, ny:{ny}")
    def _dim_check(n, W, H):
        d = H*ESC._coefficient_rectangular(s, n[0], n[1], approx=False)
        dim_x, dim_y = d*(n[0]-1), d*(n[1]-1)
        return dim_x> W[0], dim_y >W[1]

    def get_nmax(s, W, H, shape= "L", thershold = 0.3, exceed=True):
        try:
            W = list(W)
        except:
            W = [W]
        if len(W) >2:
            W = W[:2]

        if len(W) == 2 and shape=="R": #Rectangular
            # Calculate Rectagle search routine
            Wx, Wy = W
            switch = False
            m = Wy/Wx
            d = 0
            
            # initial point
            if s > 30:
                approx_d = H * ESC._coefficient_rectangular(s, N=5, M=5, approx=True)
                nx_if, ny_if = Wx/approx_d, Wy/approx_d
                n_i = [Utils.half_ceil(nx_if), Utils.half_ceil(ny_if)]
            else:
                #nx_i, m1, m2, thershold_condition = ESC._linear_nmax(s, W, H, thershold)
                #ny_i = Utils.half_ceil(m*nx_i)
                #n_i = [nx_i, ny_i]
                if m >=1:
                    n_i=[1, Utils.half_ceil(m)]
                else:
                    x = Utils.half_ceil(1/m)
                    n_i = [x if x >=1 else 1, 1]
                
            
            p_range = (5, 20)
            
            #print(f"Inital point: {n_i}")
            
            # Search
            # Line algorithm: Get adjacent points of the given ratio
            # Measuring distances from center of (Wx, Wy) plane.
            # Pick closest point to the center on (Wx, Wy) plane.
            points_generated = points(n_i, line_param=[m, d], p_range=p_range, allow_cross_thick=True)
            search_points = [p for p in points_generated if p[0]>=1 and p[1]>=1]

            # print(search_points)
            areas = [ESC._area(s, *point, shape = "R") for point in search_points]
            measures = [math.sqrt( (lx*H - Wx)**2+ (ly*H-Wy)**2) for lx, ly in areas]
            
            # print(measures)

            v, i = min([(val, index) for index, val in enumerate(measures)])
            nx, ny = search_points[i]
            esc_coef = ESC._coefficient_rectangular(s, N=nx, M=ny, approx=True)
            thershold_condition = True
            if math.fabs((nx-1)*esc_coef - Wx) > esc_coef*thershold or math.fabs((ny-1)*esc_coef - Wx) > esc_coef*thershold:
                thershold_condition = False
            N = (nx * ny, nx, ny, thershold_condition)

        elif shape=="L":
            Wx = [W[0]]
            n, m, l, ther_1 = ESC._linear_nmax(s, Wx, H, thershold)
            if len(W) ==2:
                n2, m2, l2, ther_2 = ESC._linear_nmax(s, [W[1]], H, thershold)
                N = (n*n2, n, n2, ther_1 != ther_2)
            else:
                N = (n, m, l, ther_1)

        return N
    
    def array(s, N, M=1, shape="L", approx=False):
        dx, dy = ESC.coefficient(s, N, M, shape, approx)
        dy = 0 if dy is None else dy
        xarr = dx* PositionArray.csym_index(N)
        yarr = dy* PositionArray.csym_index(M)
        return PositionArray.from_arrays(xarr, yarr)

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
            L = OP.D(x/H, W/H, s)
            try:
                sol = op.root_scalar(lambda xc: OP.D(xc, W/H, s) + L, bracket=[0,W/(2*H)], method="brentq")
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

    #--------------------------------------------------------
    def fill_rq(s, H, xm, xe, status = 0): #Done, status: 0 fill Q region including P, 1: only Q
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
        
        xe = OP.xe(s, Wx, H)
        xm = OP.xm(s, Wx, H, xe)
        ye = OP.xe(s, Wx, H)
        ym = OP.xm(s, Wx, H, ye)

        if (xe - xarr.max()) >= 0.4*dx:
            xarr[0] = xarr[N-1] = xe
        if (ye - yarr.max()) >= 0.4*dy:
            yarr[0] = yarr[M-1] = ye
        
        x_ex = OP._get_r_points(xarr, dx, s, H, Wx, xe, xm)
        y_ex = OP._get_r_points(yarr, dy, s, H, Wy, ye, ym)

        print(f"x_ex:{x_ex}")
        print(f"xarr:{xarr}")
        print(f"dx:{dx}")
        print(f"y_ex:{y_ex}")
        print(f"yarr:{yarr}")
        print(f"dy:{dy}")
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
            pass

    

class Utils: # basic Utils and functions including mathematical routines
    def d2(x, y):
        return x**2 + y**2
    def plane_meshgrid(x_range, y_range, dim):
        xline = np.linspace(x_range[0], x_range[1], dim[0])
        yline = np.linspace(y_range[0], y_range[1], dim[1])
        return np.meshgrid(xline, yline)

    def evaluation_factors(arr):
        stdE = (arr.max() - arr.min())/arr.mean() # |E_max -E_min|/E_mean
        rmse = (arr.std())/arr.mean() # std/mean
        return stdE, rmse
    def half_ceil(x):
        rho = x - math.floor(x)
        if rho >= 0.5:
            result = math.ceil(x)
        else:
            result = math.floor(x)
        return result
    def near_integers(x):
        low = math.floor(x)
        return low, low+1 
    
    def print_mesh_point(arr):
        for line in arr:
            for element in line:
                print(f"({element[0]:.2}, {element[1]:.2})", end = "")
            print(":\n")
    def plot_xy_mesh(points, scale, dim):
        scale_x, scale_y = scale
        xi, xf, yi, yf = points
        x_r = (scale_x * xi, scale_x * xf)
        y_r = (scale_y * yi, scale_y * yf)
        return Utils.plane_meshgrid(x_r, y_r, dim), [*x_r, *y_r] 
    def data_ceiling(data, n):
        if n <= 1:
            raise ValueError("n must be greater than 1 and integer.")
        d_i = data.min()
        dn = (data.max()-data.min())/n
        dn_add = (data.max()-data.min())/(n-1)
        for i in range(0, n):
            di1 = i*dn + d_i
            di2 = di1+dn
            data = np.where((data >= di1) & (data < di2 ), ((i))*dn_add + d_i, data)
        return data
    # Fix rectangle
    def rectangle_line_points(dx, dy, Wx=None, Wy=None, wx=1, wy=1):
        Wx = 0 if Wx is None else Wx
        Wy = 0 if Wy is None else Wy 
        return [
            [[wx*(Wx-dx)/2, wx*(Wx+dx)/2],[wy*(Wy-dy)/2, wy*(Wy-dy)/2]],
            [[wx*(Wx-dx)/2, wx*(Wx+dx)/2],[wy*(Wy+dy)/2, wy*(Wy+dy)/2]],
            [[wx*(Wx-dx)/2, wx*(Wx-dx)/2],[wy*(Wy-dy)/2, wy*(Wy+dy)/2]],
            [[wx*(Wx+dx)/2, wx*(Wx+dx)/2],[wy*(Wy-dy)/2, wy*(Wy+dy)/2]]
        ]

class UtilsAlgorithm:
#=============================================================================================================================================
    def transformMatrix(n):
        Fd = np.array([[1 if i == j else (-1 if i-j == 1 else 0) for j in range(0,n)] for i in range(0,n)])
        inFd =  np.array([[1 if i >= j else (0) for j in range(0,n)] for i in range(0,n)])

        return Fd, inFd
    def loc_to_diff(x, n):
        fd, infd = UtilsAlgorithm.transformMatrix(n)
        d= fd.dot(x)
        w= -2*d[0]
        d0 = d[1:]
        m = math.floor(n/2) if n%2 == 0 else math.floor((n-1)/2)
        #print(m)
        return d0[0:m], w, m
    def diff_to_loc(diff, n, w):
        fd , infd = UtilsAlgorithm.transformMatrix(n)

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
        fd, infd = UtilsAlgorithm.transformMatrix(n)
    
        if n%2 ==0:
            d0 = np.append(d[0:m], np.flip(d[0:m-1]))
            pass
        else:
            d0 = np.append(d,np.flip(d))
    
        d0 = np.insert(d0, 0, -W/2)
        xi = infd.dot(d0)
    
        location = np.array([[[x, 0] for x in xi]])
    
        return (ydata - Utils.gauss_distribution(xdata, np.array([[[0]]]), location,h))[0][0]
