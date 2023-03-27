# Python mild imeplementation of Chebfun library 2d approximation.
# Original project: https://www.chebfun.org/
# Algorithm details: Townsend, A., & Trefethen, L. N. (2013). An extension of Chebfun to two dimensions. SIAM Journal on Scientific Computing, 35(6), C495-C518.
from typing import Callable, Tuple, Union
from numbers import Number
from collections.abc import Iterable


import numpy as np
from numpy.polynomial import Chebyshev
from numpy.typing import NDArray

SingleFunction = Callable[[Number], Number]
Callable_check = np.frompyfunc(callable, 1, 1)

# Basic utils
def cheby_range_transform(x, a, b): # [-1, 1] -> [a, b]
    return ((b-a)/2)*x + (b+a)/2
def cheby_range_inv_transform(x, a, b): # [a, b] -> [-1, 1]
    return (2/(b-a))*x - (b+a)/(b-a)
def cheby_ext_grid(a, b, n):
    x_i = np.cos((np.pi/n-1)*np.arange(n))
    return cheby_range_transform(x_i, a, b)
def cheby_root_grid(a, b, n):
    x_i = np.cos((np.pi/n)*(0.5+np.arange(n)))
    return cheby_range_transform(x_i, a, b)


def _cheby_stage_1_params(i):
    i +=1
    return int(2**(i+2)+1), int(2**(i)+1)
def get_xy_cheby_decompose(
        f:Callable, 
        point:Tuple[Number, Number], 
        domain:Tuple[Number, Number, Number, Number], 
        cheby_deg: int):
    a, b, c, d = domain
    x, y = point
    # v(x)
    vx_list = cheby_root_grid(a, c, cheby_deg+1)
    vy_list = np.full_like(vx_list, y)
    vz_list = f(vx_list, vy_list)

    v_i = Chebyshev.fit(vx_list, vz_list, deg=cheby_deg, domain=(a, c))

    # u(y)
    uy_list = cheby_root_grid(b, d, cheby_deg+1)
    ux_list = np.full_like(uy_list, x)
    uz_list = f(ux_list, uy_list)

    u_i = Chebyshev.fit(uy_list, uz_list, deg=cheby_deg, domain=(b, d))
    return v_i, u_i


class RankApprox2dim:
    def __init__(self, 
                 fx_list:Union[SingleFunction, Iterable[SingleFunction]], 
                 fy_list:Union[SingleFunction, Iterable[SingleFunction]], 
                 weights:Union[Number, Iterable[Number]], 
                 w_domain = float,
                 domain:Union[None, Tuple[Number, Number, Number, Number]]=None, 
                 **kwargs):
        self._fx_list = np.array(fx_list) if not isinstance(fx_list, np.ndarray) else fx_list
        self._fy_list = np.array(fy_list) if not isinstance(fy_list, np.ndarray) else fy_list
        self._weights = np.array(weights, dtype=w_domain) if not isinstance(weights, np.ndarray) else weights

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
    def __call__(self, x:Union[Number, NDArray], y:Union[Number, NDArray], dtype:Union[None, Tuple[Union[object, type], Union[object, type]]]=None):
        if dtype is None:
            dtype = (float, float)
        if not isinstance(x, np.ndarray):
            f_x_v = np.fromiter((fi(x) for fi in self._fx_list), dtype=np.dtype(dtype[0]))
            f_y_v = np.fromiter((fi(y) for fi in self._fy_list), dtype=np.dtype(dtype[1]))
            return (self._weights * f_x_v * (f_y_v)).sum()
        # x, y vector
        # mesh = 2 dim
        shape = x.shape
        if len(shape) != 1: 
            x = x.reshape(-1)
            y = y.reshape(-1)
        if len(x) != len(y):
            raise ValueError("x, y values must have same dimension.")
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
    #-----------------------------------------------------------------------------
    @classmethod
    def from_function_approx(cls, 
                             f:Callable, 
                             domain:Tuple[Number, Number, Number, Number], 
                             tol_err = 5E-15,
                             initial_step=0,
                             max_step = 5):
        if initial_step <0:
            initial_step = 0
        if max_step < initial_step:
            raise ValueError(f"max_step < initial_step: {max_step} < {initial_step}, impossible.")
        
        x_i, x_f =domain[0], domain[2]
        y_i, y_f = domain[1], domain[3]
        
        err = 1.0
        for j in range(initial_step, max_step):
            n, max_step_stage_1 = _cheby_stage_1_params(j)

            cheby_deg= (n-1)

            arr_x = cheby_root_grid(x_i, x_f, n)
            arr_y = cheby_root_grid(y_i, y_f, n)
            arr_x.sort()
            arr_y.sort()

            p_x, p_y = np.meshgrid(arr_x, arr_y)
            _, c_l = p_x.shape

            e_k_approx = cls(lambda x: x, lambda y: y, 0)
            f_k = cls(lambda x: x, lambda y: y, 0)
            def e_k(x, y):
                return f(x, y) + e_k_approx(x, y)
            err = 1.0
            
            for k in range(0, max_step_stage_1):
                e_k_val = e_k(p_x, p_y)
                err = np.abs(e_k_val).max() 
                if err < tol_err: 
                    break
                max_index = np.argmax(e_k_val)
                r = int(max_index/c_l)
                c = max_index%c_l
                x_k, y_k = p_x[r, c], p_y[r, c]

                u_k, v_k = get_xy_decompose(e_k, (x_k, y_k ), domain, cheby_deg)
                d_k  = 1/ e_k(x_k, y_k)

                e_k_approx.rank_up(u_k, v_k, -d_k)
                f_k.rank_up(       u_k, v_k, d_k)
            
            if err < tol_err: 
                    break

        return f_k, err
    #-----------------------------------------------------------------------------
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
            self._fx_list = np.delete(self._fx_list[: -a], np.s_[-1])
            self._fy_list = np.delete(self._fy_list[: -a], np.s_[-1])
            self._weights = np.delete(self._weights[:-a] , np.s_[-1])
            pass
        elif b is None:
            # del last 'a' length elements
            self._fx_list = np.delete(self._fx_list, np.s_[-a:])
            self._fy_list = np.delete(self._fy_list, np.s_[-a:])
            self._weights = np.delete(self._weights, np.s_[-a:])
        else:
            # del from a to b elements
            self._fx_list = np.delete(self._fx_list, np.s_[a:b])
            self._fy_list = np.delete(self._fy_list, np.s_[a:b])
            self._weights = np.delete(self._weights, np.s_[a:b])
    #--------------------------------------------------------------------------------
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
