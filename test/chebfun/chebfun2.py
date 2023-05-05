from typing import Iterable, Union, Callable, Tuple
from numbers import Number
from numpy.typing import NDArray
from collections.abc import Iterable as iterable_instance

import math
import numpy as np
from numpy.polynomial.chebyshev import Chebyshev
from convolution import convolve2d

FLOAT_EPS = np.finfo(float).eps # 2.22...e-16

Number_Full= Union[Number, np.inf, np.NINF] # Common real number, infty, -infty
SingleFunction = Callable[[Number], Number]

ChebFuncPref = {
    "domain": [-1, 1],
    "splitting": 0,
    "splitPrefs" : {"splitLength": 160, "splitMaxLength":6000},
    "blowup": 0,
    "blowupPrefs":{
        "exponentTol": 1.100000e-11,
        "maxPoleorder": 20,
        "defaultSingType": 'sing'
        },
    "enableDeltaFunctions": 1,
    "deltaPrefs": {
        "deltaTol":     1.000000e-09,
        "proximityTol": 1.000000e-11
        },
    "cheb2Prefs": {
        "chebfun2eps":              2.220446e-16,
        "maxRank":                  513,
        "sampleTest":               1
        },
    "cheb3Prefs": {
        "chebfun3eps":              2.220446e-16,
        "maxRank:":                  128,
        "sampleTest":               1,
        "constructor":              "chebfun3f"
        },
    "tech":                         "@chebtech2",
    "techPrefs": {
        "chebfuneps":              2.220446049250313e-16,
        "minSamples":               17,
        "maxLength":                65537,
        "fixedLength":              np.NaN,
        "extrapolate":              0,
        "sampleTest":               1,
        "refinementFunction":       'nested',
        "happinessCheck":           'standard',
        "useTurbo":                 0
        }
}



def cheby_range_transform(x, a, b): # [-1, 1] -> [a, b]
    return ((b-a)/2)*x + (b+a)/2
def cheby_range_inv_transform(x, a, b): # [a, b] -> [-1, 1]
    return (2/(b-a))*x - (b+a)/(b-a)
def cheby_ext_grid_1(a, b, n):
    x_i = np.cos((np.pi/n-1)*np.arange(n))
    return cheby_range_transform(x_i, a, b)
def cheby_root_grid_1(a, b, n): # Same with "chebypts(n, (a,b), 1)""
    x_i = np.cos((np.pi/n)*(0.5+np.arange(n)))
    return cheby_range_transform(x_i, a, b)

def chebpts(n, domain=[-1, 1], pol_type = 2):
    if pol_type ==2:
        x = np.cos(np.pi/(n-1)*np.arange(n)) if n>1 else np.ones(1)
    elif pol_type ==1:
        x = np.cos((np.pi/n)*(0.5+np.arange(n)))
    return cheby_range_transform(x, domain[0], domain[1])

def trigpts(n, dom):
    x_i, x_f = dom
    size =(x_f - x_i )/2
    center = (x_f+ x_i)/2
    x = np.linspace(-np.pi, np.pi, n+1)/np.pi
    x_center = (x[0] + x[-1])/2
    return size*(x-x_center) + center


def _ge_pivots(vals:NDArray, tol_abs:Number, factor:Number): #= completeACA

        [nx, ny] = vals.shape
        width = min(nx, ny)
        pivots_value = [0]
        pivots_elements = [[0,0]]
        ifail = 1

        z_rows = 0
        ind = np.argmax(vals)
        row, col = divmod(ny, ind) # ind2suv
        infNorm = np.abs(vals[row, col])

        # Diagonal bias 
        if nx == ny and np.max(np.abs(vals.diagonal())) - infNorm > - abs_tol:
            ind = np.argmax(np.abs(vals.diagonal()))
            row = ind
            col = ind
            infNorm = np.abs(vals[row,col])

        scl = infNorm

        if np.abs(scl - 0) < FLOAT_EPS:
            pivot_value = 0.
            rows = np.zeros()
            cols = np.zeros()
            ifail = 0
        else:
            rows = np.zeros((1, vals.shape[1]))
            cols = np.zeros((vals.shape[0], 1))
            

        while (infNorm > tol_abs and z_rows < width/factor ):
            row_vec = vals[:, col].reshape(-1, 1)
            col_vec = vals[row ,:].reshape(1, -1)
            pivot_val = vals[row, col]
            dvals = (convolve2d(row_vec, col_vec))/pivot_val
            vals = vals - dvals
            
            z_row += 1
            row ,col = divmod(np.argmax(np.abs(vals)), vals.shape[1])

            pivot_value.append(pivot_val)
            pivots_elements.append([row, col])

            if nx == ny and np.max(np.abs(vals.diagonal())) - infNorm > - tol_abs:
                ind = np.argmax(np.abs(vals.diagonal()))
                row = ind
                col = ind
                infNorm = np.abs(vals[row,col])

        ifail = 0 if infNorm <= tol_abs else ifail
        ifail = 1 if z_rows >= width/factor else ifail 

        return np.array(pivot_value), np.array(pivots_elements), rows, cols, ifail
def _get_tol(X, Y, vals, domain, pseudo_level):
        #Done
        nx, ny = vals.shape
        grid = max(nx, ny)
        dfdx = 0
        dfdy = 0

        if nx > 1 and ny >1 :
            dfdy = np.diff(vals[:, :], axis=1) / np.diff(Y[:-1, :], axis=1)
            dfdx = np.diff(vals[:, :], axis=0) / np.diff(X[:, :-1], axis=0)
        elif nx >1 and ny ==1:
            dfdx = np.diff(vals, axis=1) / np.diff(X, axis=1)
        elif nx == 1 and ny>1:
            dfdy = np.diff(vals, axis=0) / np.diff(Y, axis=0)
        
        jac_norm = np.max(np.maximum(np.abs(dfdx), np.abs(dfdy)))
        vscale = np.max(np.abs(vals))
        rel_tol = grid**(2/3) * pseudo_level
        abs_tol = np.max(np.abs(domain)) * np.max(jac_norm, vscale) * rel_tol
        return rel_tol, abs_tol
def _point_mesh(nx:int, ny:int, dom = [-1, 1, -1, 1], fx_type=2, fy_type=2):
    # type: 0: tri, 1: T of 1st kind, 2: U, of 2nd kind
    x_i, x_f, y_i, y_f = dom
    x_arr = chebpts(nx, (x_i, x_f), pol_type=fx_type) if fx_type != 0 else trigpts(nx, (x_i, x_f)) 
    y_arr = chebpts(ny, (y_i,y_f), pol_type=fy_type) if fx_type != 0 else trigpts(nx, (y_i, y_f)) 
    return np.meshgrid(x_arr, y_arr)

def _grid_refine(n:int, f_type:int):
    if f_type == 0: # tri
        grid =  2**(math.floor(math.log2(n)+1))
        nesting = np.arange(0, grid, 2) # matlab "1: 2: grid"
    elif f_type == 1:
        grid = 3*n
        nesting = np.arange(1, grid, 3)
    elif f_type == 2:
        grid = 2**(math.floor(math.log2(n)+1)) +1
        nesting = np.arange(0, grid, 2)

    return grid, nesting



class ChebFun2: # Only support smooth function
    factor = 4 # ratio between size of matrix and number of pivots.
    def __init__(self, 
                 func_x=Union[None, Tuple[Callable]], 
                 func_y=Union[None, Tuple[Callable]], 
                 weights=None,
                 domain = [-1, 1, -1, 1]):
        self.func_x = [func_x] if func_x is not None else [] # col vector
        self.func_y = [func_y] if func_y is not None else [] # row vector
        
        n = len(self.func_x)
        if weights is None:
            self.weights = np.array([], dtype=float)
        elif isinstance(weights, (iterable_instance, Number)):
            weights = np.array(weights)
            s = weights.shape
            ns = len(s)
            if ns >2 :
                weights = weights.reshape((n,n))
            if ns ==1:
                self.weights = np.zeros((n,n))
                np.fill_diagonal(self.weights, weights)
            elif ns ==2:
                self.weights = weights 
    @property
    def rank(self):
        return np.max(self.weights.shape)
    def __str__(self):
        return "Chebfun2 object, domain:{self.domain}, rank:{self.rank}, eps:{self.eps}"
    def __call__(self, x, y): # Support numpy mesh array
        if self.func_x is None:
            return 0 if not isinstance(x, np.ndarray) else np.zeros(x.shape)
        if not isinstance(x, np.ndarray):
            f_x_v = np.fromiter((f(x) for f in self.func_x), dtype=np.dtype(float))
            f_y_v = np.fromiter((f(y) for f in self.func_y), dtype=np.dtype(float))
            return np.einsum("i,ij,j", f_y_v, self.weights, f_x_v)
        # 2 dim mesh
        shape = x.shape
        if len(shape) != 1: 
            x = x.reshape(-1)
            y = y.reshape(-1)
        if len(x) != len(y):
            raise ValueError("x, y values must have same dimension.")
        fx_v = np.vstack(f(x) for f in self.func_x)
        fy_v = np.vstack(f(y) for f in self.func_y)
        try:
            result = np.einsum("ji,jk,ki->i",fy_v, self.weights, fx_v).reshape(shape)
        except:
            print(fy_v.shape, self.weights.shape, fx_v.shape)
            raise ValueError("See above message")
        return result
    def __add__(self, outer):
        pass
    def __sub__(self, outer):
        pass
    def __mul__(self, outer):
        pass
    def diff(self, axis=None):
        pass
    def integral(self, axis=1):
        pass
    @classmethod
    def constructor(cls, 
                 func:Callable, 
                 domain:Tuple[float, ...], # infty range can be passed in 
                 istri=False,
                 isequi = False,
                 samples:Tuple[int, ...] = [17, 17],
                 max_rank =None,
                 float_eps = None #=pseudo_level
                 ):
        if callable(func):
            if isinstance(func, ChebFun2): # copy the objects
                pass
            pass
        elif isinstance(func, iterable_instance): # Assum all elements are callable
            pass

        if float_eps is None:
            float_eps = FLOAT_EPS

        if len(samples) > 4:
            samples = samples[:4]
        
        if len(samples) == 3:
            samples.append(samples[-1])
        elif len(samples) == 2:
            samples.append()
        
        if len(samples) == 4 :
            min_sample = samples[:2]
            max_sample = samples[2:4]

        prefx=  {}
        prefy= {}
            
        resolved = False
        failed = False
        min_sample:Tuple[int, int] =  2**(np.log2()) + (0 if istri else 1) 
        #--------------------------------
        while(not resolved and not failed):
            grid = min_sample
            X, Y = _point_mesh(*grid, domain, prefx, prefy)
            vals = func(X, Y)

            vscale = np.max(np.abs(vals))
            if np.isinf(vscale):
                raise RuntimeError(f"The given function has some singular values in the given domain:{domain}")
            elif np.isnan(vals):
                raise RuntimeError(f"The given function failed to evaluate function value in the given fomain {domain}")
            
            
            tol_rel, tol_abs = _get_tol(X, y, vals, domain, float_eps)
            
            # Stage 1: Pivot searching
            [  pivot_values, pivot_position, 
                func_y, func_x, ifail
            ] = _ge_pivots(vals, tol_abs, cls.factor)
            strike = 1
        
            while(
                resolved and  (grid <= cls.factor*(max_rank -1)+1) and strike < 3
                ):
                grid = _grid_refine(grid, prefx, prefy)
                
                X, Y = _point_mesh(*grid, domain, prefx, prefy)
                vals = func(X,Y)
                vscale = np.max(np.abs(vals))

                tol_rel, tol_abs=  _get_tol(X, Y, vals, domain, float_eps)

                [  pivot_values, pivot_position, 
                    func_y, func_x, ifail
                ] = _ge_pivots(vals, tol_abs, cls.factor)
                if np.abs(pivot_values) < 1E4*vscale*tol_rel:
                    strike +=1
            
            if np.any(grid > cls.factor * (max_rank-1) + 1):
                raise RuntimeWarning("Not a low-rank function")
                failed = True
            
            # Pivot slices resolve checking


            # Stage 2

            



            rel_tol, abs_tol = self._get_tol(X, Y, Vals, domain, pseudo_level=float_eps)

        resolved # = isHappy

        # Stage 2: Resolve column and row slices
        pass
    
    @classmethod
    def from_numeric_data(self, data:NDArray, domain):
        pass


#---------------------------------------------------

# Old =================
class RankApprox2dim:
    def __init__(self, 
                 fx_list:Union[SingleFunction, Iterable[SingleFunction]], 
                 fy_list:Union[SingleFunction, Iterable[SingleFunction]], 
                 weights:Union[Number, Iterable[Number]], 
                 weight_dtype = float,
                 domain:Union[None, Tuple[Number, Number, Number, Number]]=None, 
                 **kwargs):
        self._fx_list = np.array(fx_list) if not isinstance(fx_list, np.ndarray) else fx_list
        self._fy_list = np.array(fy_list) if not isinstance(fy_list, np.ndarray) else fy_list
        self._weights = np.array(weights, dtype=weight_dtype) if not isinstance(weights, np.ndarray) else weights

        if (self._fx_list.shape) != 1:
            self._fx_list = self._fx_list.reshape(-1) 
        if (self._fy_list.shape) != 1:
            self._fy_list = self._fy_list.reshape(-1) 
        if (self._weights.shape) != 1:
            self._weights = self._weights.reshape(-1) 

        l_a = len(self._fx_list)
        l_b = len(self._fy_list)
        l_c = len(self._weights)

        if l_a != l_b or l_b != l_c or l_a != l_c:
            raise ValueError("Dimensions are not same, {l_a}, {l_b}, {l_c}.")

        self.__check_types(self._fx_list, self._fy_list, self._weights)
        self.domain = domain
        self.add_info = kwargs
    def __add__(self, other):
        x_i1, y_i1, x_f1, y_f1 = self.domain
        x_i2, y_i2, x_f2, y_f2 = other.domain

        x_i = x_i1 if x_i1 > x_i2 else x_i2
        x_f = x_f1 if x_f1 < x_f2 else x_f2
        y_i = y_i1 if y_i1 > y_i2 else y_i2
        y_f = y_f1 if y_f1 < y_f2 else y_f2
        domain = [x_i, y_i, x_f, y_f]

        fx_list = np.concatenate([self.fx_list, other.fx_list])
        fy_list = np.concatenate([self.fy_list, other.fy_list])
        weights = np.concatenate([self.weights, other.weights])
        return RankApprox2dim(fx_list, fy_list, weights, domain= domain)
    def __sub__(self, other):
        x_i1, y_i1, x_f1, y_f1 = self.domain
        x_i2, y_i2, x_f2, y_f2 = other.domain

        x_i = x_i1 if x_i1 > x_i2 else x_i2
        x_f = x_f1 if x_f1 < x_f2 else x_f2
        y_i = y_i1 if y_i1 > y_i2 else y_i2
        y_f = y_f1 if y_f1 < y_f2 else y_f2
        domain = [x_i, y_i, x_f, y_f]

        fx_list = np.concatenate([self.fx_list, other.fx_list])
        fy_list = np.concatenate([self.fy_list, other.fy_list])
        weights = np.concatenate([self.weights, -other.weights])
        return RankApprox2dim(fx_list, fy_list, weights, domain= domain)
    def __mul__(self, other:Number): #Scalar
        weights= other*self._weights
        return RankApprox2dim(self._fx_list, self._fy_list, weights, domain = self.domain)
    def __truediv__(self, other:Number):
        weights= self._weights/other
        return RankApprox2dim(self._fx_list, self._fy_list, weights, domain = self.domain)
    def diff(self, axis = 2):
        pass
    def integral(self, axis =2):
        pass
    def __call__(self, x:Union[Number, NDArray], y:Union[Number, NDArray], dtype:Union[None, Tuple[Union[object, type], Union[object, type]]]=None):
        if dtype is None:
            dtype = (float, float)
        if not isinstance(x, np.ndarray):
            f_x_v = np.fromiter((fi(x) for fi in self._fx_list), dtype=np.dtype(dtype[0]))
            f_y_v = np.fromiter((fi(y) for fi in self._fy_list), dtype=np.dtype(dtype[1]))
            return (self._weights * f_x_v * (f_y_v)).sum()
        # x, y vector
        # mesh = 2 dim
        shape = x.shape if x.size > y.size else y.shape
        if len(shape) != 1: 
            x = x.reshape(-1)
            y = y.reshape(-1)
        if x.size != y.size:
            raise ValueError(f"x, y values must have same dimension., {x.size} {y.size}")
        weight = np.tile(self._weights, [len(x), 1]).transpose()

        fx_v = np.vstack(f(x) for f in self._fx_list)
        fy_v = np.vstack(f(y) for f in self._fy_list)
        return (weight * fx_v * fy_v).sum(axis=0).reshape(shape)
    def __check_types(self, fx_list:NDArray, fy_list:NDArray, weights:NDArray):
        if not Callable_check(fx_list).all():
            raise TypeError("fx list must consist of callable objects.")
        if not Callable_check(fy_list).all():
            raise TypeError("fy list must consist of callable objects.")
        if  not np.issubdtype(weights.dtype, np.number):
            raise TypeError("Weight must consist of number type.") 
    @classmethod
    def from_function_approx(cls, 
                             f:Callable, 
                             domain:Tuple[Number, Number, Number, Number], 
                             tol_err = 5E-15,
                             initial_step=0,
                             max_step = 6,
                             include_boundary= False,
                             hold_pivot = False,
                             show_process=False):
        if initial_step <0:
            initial_step = 0
        if max_step < initial_step:
            raise ValueError(f"max_step < initial_step: {max_step} < {initial_step}, impossible.")
        
        x_i, x_f =domain[0], domain[2]
        y_i, y_f = domain[1], domain[3]
        
        err = 1.0
        
        # Stage 1
        for j in range(initial_step, max_step):
            n, max_step_stage_1 = Cheby_stage_params(j)

            cheby_deg= (n-1) if not include_boundary else n+1

            arr_x = cheby_root_grid(x_i, x_f, n)
            arr_y = cheby_root_grid(y_i, y_f, n)
            if include_boundary:
                arr_x = np.insert(arr_x, -1, [x_i, x_f])
                arr_y = np.insert(arr_y, -1, [y_i, y_f])
            arr_x.sort()
            arr_y.sort()

            p_x, p_y = np.meshgrid(arr_x, arr_y)
            _, c_l = p_x.shape

            #e_k_approx = cls(lambda x: x, lambda y: y, 0)
            e_k_approx = GeApprox()
            #f_k = cls(lambda x: x, lambda y: y, 0)
            def e_k(x, y):
                return f(x, y) + e_k_approx(x, y)
            err = 1.0
            pivots = []
            
            for k in range(0, max_step_stage_1):
                e_k_val = np.abs(e_k(p_x, p_y))
                err = e_k_val.max() 
                if show_process:
                    print(100*" ", end="\r")
                    print(f"Step: {j}, Rank:{e_k_approx.rank}, Err:{err}", end="\r")
                if err < tol_err: 
                    break
                max_index = np.argmax(e_k_val)
                r = int(max_index/c_l)
                c = max_index%c_l
                x_k, y_k = p_x[r, c], p_y[r, c]
                pivots.append([x_k, y_k])

                u_k, v_k = get_xy_decompose(e_k, (x_k, y_k ), domain, cheby_deg)
                
                d_k  = 1/ e_k(x_k, y_k)

                #e_k_approx.rank_up(u_k, v_k, -d_k)
                e_k_approx.add_term(u_k, v_k, -d_k)
                #f_k.rank_up(       u_k, v_k, d_k)
            
            if err < tol_err: 
                    if show_process:
                        print("")
                    break
        
        #f_k.add_infos(err = err, pivots = np.array(pivots)) if hold_pivot else f_k.add_infos(err = err)
        e_k_approx.weights = - e_k_approx.weights
        return e_k_approx
    
    def rank_up(self, 
                        fx_list:Union[SingleFunction, Iterable[SingleFunction]], 
                        fy_list:Union[SingleFunction, Iterable[SingleFunction]], 
                        weights:Union[Number, Iterable[Number]]):
        
        fx_list = list(fx_list) if isinstance(fx_list, Iterable) and hasattr(fx_list, "__getitem__") else [fx_list]
        fy_list = list(fy_list) if isinstance(fy_list, Iterable) and hasattr(fy_list, "__getitem__") else [fy_list]
        weights = list(weights) if isinstance(weights, Iterable) else [weights]

        l_a = len(fx_list) 
        l_b = len(fy_list)
        l_c = len(weights)

        if l_a != l_b or l_b != l_c or l_a != l_c:
            raise ValueError(f"Dimensions are not same, {l_a}, {l_b}, {l_c}.")
        
        self.__check_types(fx_list, fy_list, np.array(weights))

        self._fx_list = np.concatenate([self._fx_list, fx_list])
        self._fy_list = np.concatenate([self._fy_list, fy_list])
        self._weights = np.concatenate([self._weights, weights])
    def rank_down(self, a:None, b:None):
        if a is None and b is None:
            # del last element
            pass
        elif b is None:
            # del last 'a' length elements
            pass
        else:
            # del from a to b elements
            pass

    #--------------------------------------------------------------------------------
    def add_infos(self, **kwargs):
        for key in kwargs.keys():
            self.add_info[key] = kwargs[key]
    
    def elements(self, a:int, b:Tuple[None, int] = None):
        if b is None: 
            if a > self.rank-1:
                raise IndexError
        else:
            if a > b :
                a, b = b,a
            if b > self.rank-1:
                raise IndexError
        fx_e = self.fx_list[:a] if b is None else self.fx_list[a:b]
        fy_e = self.fy_list[:a] if b is None else self.fy_list[a:b]
        weight_e = self.weights[:a] if b is None else self.weights[a:b]
        return np.vstack([fx_e, fy_e, weight_e]).transpose()
    def element(self, i):
        if i > self.rank-1:
            raise IndexError
        return np.array([self.fx_list[i], self.fy_list[i], self.weights[i]])
    def partial(self, a:int, b:int = None):
        if b is not None:
            fx, fy, weight = self.elements(a, b).transpose()
        else:
            fx, fy ,weight = self.element(a)
        return RankApprox2dim(fx , fy, weight, domain=self.domain, **self.add_info)
    @property
    def fx_list(self):
        return self._fx_list
    @fx_list.setter
    def fx_list(self, value):
        value = np.array(value) if isinstance( value, np.ndarray) else value
        if len(value.shape) !=1:
            value = value.reshape(-1)
        if value.size != self._fx_list.size:
            raise ValueError(f"Length of the given function array is not same with rank {self._fx_list.size}")
        for fx in value:
            if isinstance(fx, Callable):
                continue
            raise TypeError("The given array must consist of callable objects.")
        self._fx_list = value
    @property
    def fy_list(self):
        return self._fy_list
    @fy_list.setter
    def fy_list(self, value):
        value = np.array(value) if isinstance( value, np.ndarray) else value
        if len(value.shape) !=1:
            value = value.reshape(-1)
        if value.size != self._fy_list.size:
            raise ValueError(f"Length of the given function array is not same with rank {self._fy_list.size}")
        for fy in value:
            if isinstance(fy, Callable):
                continue
            raise TypeError("The given array must consist of callable objects.")
        self._fy_list = value  
    @property
    def weights(self):
        return self._weights
    @weights.setter
    def weights(self, value):
        value = np.array(value) if isinstance( value, np.ndarray) else value
        if len(value.shape) !=1:
            value = value.reshape(-1)
        if value.size != self._weights.size:
            raise ValueError(f"Length of the given weights array is not same with rank {self._fy_list.size}")
        for w in value:
            if isinstance(w, Number):
                continue
            raise TypeError("The given array must consist of number objects.")
        self._weights = value 
    @property
    def rank(self):
        return self._weights.size
class GeApprox:
    def __init__(self, func_x=None, func_y=None, weights=None):
        self.func_x = [func_x] if func_x is not None else None
        self.func_y = [func_y] if func_y is not None else None
        if func_x is not None:
            n = len(self.func_x)
            self.weights = np.array(weights).reshape(n, n) if weights is not None else None
        else:
            self.weights = np.array([[weights]])    
    def __str__(self):
        return f"""sum_i sum_j w_ij u_i(x) v_j(y)
        {self.func_y}
        {self.weights}
        {self.func_x}"""   
    def __call__(self, x, y):
        if self.func_x is None:
            return 0 if not isinstance(x, np.ndarray) else np.zeros(x.shape)
        if not isinstance(x, np.ndarray):
            f_x_v = np.fromiter((f(x) for f in self.func_x), dtype=np.dtype(float))
            f_y_v = np.fromiter((f(y) for f in self.func_y), dtype=np.dtype(float))
            return np.einsum("i,ij,j", f_y_v, self.weights, f_x_v)
        # 2 dim mesh
        shape = x.shape
        if len(shape) != 1: 
            x = x.reshape(-1)
            y = y.reshape(-1)
        if len(x) != len(y):
            raise ValueError("x, y values must have same dimension.")
        fx_v = np.vstack(f(x) for f in self.func_x)
        fy_v = np.vstack(f(y) for f in self.func_y)
        try:
            result = np.einsum("ji,jk,ki->i",fy_v, self.weights, fx_v).reshape(shape)
        except:
            print(fy_v.shape, self.weights.shape, fx_v.shape)
            raise ValueError(" See above message")
        return result
    def add_term(self, func_x, func_y, weights):
        if self.func_x is None:
            # scalar case
            self.func_x = [func_x]
            self.func_y = [func_y]
            self.weights = np.array([[weights]]) if  not isinstance(weights, np.ndarray) else weights 
            return 0
        if isinstance(func_x, Callable):
            weights_pre = np.pad(self.weights, ((0,1),(0,1)), mode="constant", constant_values=(0,0))
            if not isinstance(weights, np.ndarray) and isinstance( weights, Number):
                n = len(self.func_x)
                weights_ = np.zeros((n+1, n+1))
                weights_[n,n] = weights
                weights = weights_
            if weights_pre.shape != weights.shape:
                raise ValueError("Invaild weights dimension, must be a (n+1) square matrix.")
            self.weights = weights_pre + weights
            self.func_x.append(func_x)
            self.func_y.append(func_y)
            return 0
        
        raise NotImplementedError("Not implemented")
    @property
    def rank(self):
        return len(self.func_x) if self.func_x is not None else 0
class GeFun(GeApprox):
    def __init__(self, func:Callable, func_x=None, func_y=None, weights=None):
        super().__init__(func_x, func_y, weights)
        self.func = func
        self.points = None
        self.mode = True
    def __call__(self, x,y):
        r1 = super().__call__(x,y)
        r2 = self.func(x,y) if self.mode else 0
        return r1 + r2
    def add_point(self, point, neg = False):
        x_k, y_k = point

        denominator = self(x_k, y_k) 
        denominator = - denominator if neg else denominator

        if self.points is not None:
            u_k_x_val = np.fromiter((u(x_k) for u in self.func_x), dtype=np.dtype(float))
            v_k_y_val = np.fromiter((v(y_k) for v in self.func_y), dtype=np.dtype(float))
            w_x = np.concatenate([self.weights @ u_k_x_val, np.ones(1)])
            w_y = np.concatenate([v_k_y_val @ self.weights, np.ones(1)])
            weights_W = convolve2d(w_x.reshape(1, -1).T, w_y.reshape(1, -1))/denominator
        else:
            self.points = np.array([point])
            weights_W = np.array([[1/denominator]])

        u_x = partial(self.func, y=y_k)
        v_y = partial(self.func, x_k)

        self.points = np.concatenate([self.points, [point]])
        self.add_term(u_x, v_y, weights = weights_W)
    def get_approximation(self):
        func = GeApprox()
        func.func_x = self.func_x
        func.func_y = self.func_y
        func.weights = -self.weights
        return func
   

# Chebyshev----------------------------------------------------

