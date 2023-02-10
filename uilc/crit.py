import math
from typing import Tuple, Literal, Union

from scipy.optimize import root_scalar, minimize_scalar

from uilc import PositionArray
from uilc.utils import float_eps # constants
from uilc.utils import d2, half_ceil, csym_index # methods
from uilc.utils import bresenham # module

# Expanded Sparrow Criterion
# Implementation and additional utils
# Reference: 

class ESC:
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
    

    @classmethod
    def _coef_linear(cls, s:float, N:int, approx:bool = False)->Tuple[float, float]: #_coefficient_linear
        if N == 2:
            cof = math.sqrt(4/(s+3))
        elif approx and (N >4 and s >30):
            cof = math.sqrt(3.2773/(s+4.2539))
        elif N%2 == 0 :
            sol = root_scalar(lambda D: cls._op_func_linear(D, s, N), bracket=[0,1], method = "brentq")
            cof = sol.root
        else:
            res = minimize_scalar(lambda D: cls._op_func_linear(D, s, N), bounds=(0,1), method = "bounded")
            cof = res.x
        return cof
    @classmethod
    def _coef_rectangular(cls, s:float, N:int, M:int, approx:bool=False)->Tuple[float, float]:
        if M > N:
            N, M = M, N
        if N==2 and N == M:
            cof= math.sqrt(4/(s+2))
        if approx == True and (N > 4 and M > 4 and s>30):
            cof= math.sqrt(1.2125/(s-3.349))
        else:
            try:
                sol = root_scalar(lambda D: cls._op_func_rectangularangular(D, s, N, M), bracket=[0,1],method="brentq")
                if sol.converged == False:
                    res = minimize_scalar(lambda D: cls._op_func_rectangularangular(D, s, N, M), bounds=(0,1), method = "bounded")
                    cof= res.x
                else:
                    cof= sol.root
            except:
                res = minimize_scalar(lambda D: cls._op_func_rectangularangular(D, s, N, M), bounds=(0,1), method = "bounded")
                cof= res.x
        return cof
    @classmethod
    def coefficient(
        cls,
        s:float,
        N:int, M:int=1,
        shape:Literal["L", "R"]="L",
        approx:bool=False
    )-> Tuple[float, float]:
        cls._parameter_check(s, N, M)
        N = half_ceil(N)
        M = half_ceil(M)

        if shape == "L":
            cof_x = cls._coef_linear(s, N, approx)
            cof_y = None if M == 1 else cls._coef_linear(s, M, approx)
        elif shape == "R":
            cof_x = cof_y = cls._coef_rectangular(s, N, M, approx)
        else:
            raise ValueError("\"shape\" argument must be \"L\" or \"R\" current value is {}".format(shape))
        
        return cof_x, cof_y

    @classmethod
    def area(cls, s:float, nx:int, ny:int=1, shape:Literal["L", "R"]="L"):
        cls._parameter_check(s, nx, ny)
        dx, dy = cls.coefficient(s, nx, ny)
        if shape == "L":
            return (nx-1)*dx, 0
        elif shape == "R" and ny >=1:
            return (nx -1)*dx, (ny-1)*dy
        else:
            raise ValueError(f"shape:{shape}, ny:{ny}")
    @classmethod
    def _linear_nmax(cls, 
            s:float, W:float, H:float, 
            thershold:float=0.3, permit_exceed=True
            )->Union[int, None]:
        xlim =W

        n_2 = 2
        length_2 = H*cls.area(s, n_2)[0]
        n_3 = 3
        length_3 = H*cls.area(s, n_3)[0]

        if length_2 < xlim:
            n = n_2
        elif length_3 < xlim:
            n = n_3
        else:
            n = n_3
            length = length_3
            while (length < xlim):
                n += 1
                length = H*cls.area(s, n_3)[0]
        
        # Choose a closest number to boundary
        if n > 3:
            ns = [n-2, n-1, n, n+1]
        else:
            if n == 3: #check 2, 3
                ns = [2, 3, 4, 5]
            else:
                ns = [2, 3, 4]

        ds = [ H*cls.coefficient(s, n)[0] for n in ns]
        residual =[(W - d*(n-1))/2 for d, n in zip(ds, ns)]
        # check exceeding permition
        ns = [n for n, resi in zip(ns, residual) if resi > 0 or permit_exceed]
        ds = [d for d, resi in zip(ds, residual) if resi >0 or permit_exceed]
        residual =[resi for resi in residual if resi>0 or permit_exceed]

        thersholds = [ math.fabs(resi) < thershold*d for d, resi in zip(ds, residual)]
        
        n_ther = [ n for n, thers in zip(ns, thersholds) if thers ]
        n_not_ther = [ n for n, thers in zip(ns, thersholds) if not thers ]
        
        if len(n_ther) != 0:
            N = max(n_ther)
            ther = True
        elif len(n_not_ther) != 0:
            N = n-1 if (n-1) in n_not_ther else n
            ther = False
        else:
            raise RuntimeError("Impossible")
        return (N, N, 1, ther)
    # Rectangular uniform distribution search
    @classmethod
    def _rect_point_clasification(cls, p, s, Wx, Wy):
        nx, ny = p
        area = cls.area(s, nx, ny, shape="R")
        x_b, y_b = area[0] > Wx, area[1] > Wy
        if x_b == y_b:
            if x_b and y_b: # > , >
                status = -1
            else: # <, <
                status = 1
        else:
            status = 0
        return status
    @classmethod
    def _rect_point_residual(cls, p, s, Wx, Wy):
        nx, ny = p
        area = cls.area(s, nx, ny, shape="R")
        return (Wx - area[0], Wy - area[1])

    @classmethod
    def _rectangular_nmax(
        cls, 
        s:float, 
        W:Tuple[float, float], H:float, 
        thershold:float=0.3, permit_exceed=True,
        iter_max = 1000
        ):
        # This function search point range 
        # using Bresenham line algorithm, "bresenham" module.

        Wx, Wy = W
        m = Wy/Wx
        d = 0
        Wx, Wy = Wx/H, Wy/H
        

        # initial point setting
        if s>30:
            approx_d = cls.coefficient(s, N=5, M=5, approx=True)
            p_i = [half_ceil(Wx/approx_d), half_ceil(Wy/approx_d)]
        else:
            # Using line nmax
            n_x = cls._linear_nmax(s, Wx, H, permit_exceed)
            p_i = [n_x, half_ceil(m * n_x)] 
        
        # status
        a = 1 # both larger than hor and ver values of the area
        b = 0 # larger, smaller or smaller, larger
        c = -1 # both smaller than the hor and ver values of the area       
        residual_i = cls._rect_point_residual(p_i, s, Wx, Wy)
        status_i = cls._rect_point_clasification(p_i, s, Wx, Wy)
        
        # Search
        # # 1. Get adjacent points near 'y = mx' line for given 'x'.
        # # 2. Calculate measure of (N_i, M_i) points at 1
        # # 3. If exceeding event is occured, pick current and previous point.
        point_range = [None, None]
        p = p_i
        status = status_i
        residual = residual_i
        dir = True if status_i != c else False

        iteration_max = iter_max
        iter_i = 0
        while(iter_i < iteration_max):
            # Get new points and estimate them
            p_current, p2, p3 = bresenham.next_points(p, m, dir)
            status_current = cls._rect_point_clasification(p_current, s, Wx, Wy)
            residual_current = cls._rect_point_residual(p_current, s, Wx, Wy)
            # Break, or manipulate direction
            if status_current != status:
                if status_current == b:
                    if status == a:
                        point_range[0] = p
                    else:
                        point_range[1] = p
                else: # direction change
                    if status == -status_current:
                        pi, pf = p, p_current if p_current[0] > p[0] else p_current, p 
                        point_range = [pi, pf]
                        break
                    elif status == b:
                        if status_current == a:
                            point_range[0] = p_current
                        else:
                            point_range[1] = p_current
                    
                    dir = not dir
                    continue

            else: # Check that residual is decreasing or increasing.
                if status != b: # a and c states, same if "status_current != b"
                    dfx, dfy = residual_current[0] - residual[0], residual_current[1] - residual[1]
                    if dfx*dfy >0 and status * dfx >0: 
                        # ((fx'>0 and fy'>0) or (fx'<0 and fy'<0)) 
                        # and ((a:1 and f'>0) or (c:-1 and f'<0))
                        dir = not dir
                        continue

            if point_range[0] is not None and point_range[1] is not None:
                break

            p = p_current
            status = status_current
            residual = residual_current
            
        if iter_i >= iteration_max and point_range[0] is None:
            raise RuntimeError(f"Cannot find point in {iteration_max} iteration.")

        # Get estimation range
        pi, pf = point_range 
        points_main, points_sub = bresenham.points(
            pi, line_param = [m, d], 
            p_range=[0, pf-pi], allow_corss_p = True)
        # Estimation - in progress
        if not permit_exceed:
            residuals_main =[cls._rect_point_residual(p, s , Wx, Wy), for p in points_main]
            residuals_sub = [cls._rect_point_residual(p, s , Wx, Wy), for p in points_sub]
            points_main = [p for p, resi in zip(points_main, residuals_main) if resi[0]>0 and resi[1]>0]
            points_sub  = [p for p, resi in zip(points_sub , residuals_sub)if resi[0]>0 and resi[1]>0]
        residual_norm_main = [d2(p[0], p[1]) for p in points_main]
        residual_norm_sub = [d2(p[0], p[1]) for p in points_sub]

        nx, ny = p_final
        return nx*ny, nx, ny, thers
    @classmethod
    def nmax_for_region(cls,
        s:float, 
        W:Union[float, Tuple[float, float]], H:float, 
        shape:Literal["L", "R"] = "L", 
        thershold:float = 0.3,
        permit_exceed:bool = True
        ) -> Tuple[int, int, int, bool]:
        if hasattr(W, "__iter__"):
            if len(W) >2:
                W = W[:2]
                Wx, Wy = W
        else:
            Wx = W
            Wy = None

        if shape == "L":
            n1, m1, l1, ther_1 = cls._linear_nmax(s, Wx, H, thershold, permit_exceed)
            N = (n1, m1, l1, ther_1)
            if Wy is not None:
                n2, m2, l2, ther_2 = cls._linear_nmax(s, Wx, H, thershold, permit_exceed)
                N = (m1 * m2, m1, m2, ther_1 and ther_2)
        elif shape =="R":
            if Wy is None:
                Wy = Wx
            switch = False
            if Wy > Wx:
                Wx, Wy = Wy, Wx
                switch = True
            nm, n, m, ther = cls._rectangular_nmax(s, (Wx, Wy), H, thershold, permit_exceed)
            N = (nm, m, n , ther) if switch else (nm, n, m, ther)
        else:
            raise ValueError(f"\"shape\" must be \"L\" or \"R\" current:{shape}")
    
    @classmethod
    def array(cls, s:float, N:int, M=1, shape="L", approx=False):
        cls._parameter_check(s, N, M)

        dx, dy = cls.coefficient(s, N, M, shape, approx)
        
        xarr = dx*csym_index(N)
        
class BC:
`   `