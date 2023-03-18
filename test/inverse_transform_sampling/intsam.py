import math

from typing import Callable, Tuple, Union
from numpy.typing import NDArray, ArrayLike
from numbers import Number
from collections.abc import Iterable

import numpy as np

from numpy.polynomial.chebyshev import Chebyshev
from scipy.optimize import bisect

EPS = 10E6*np.finfo(float).eps

# Discrete version
def pmf2cdf(prob_mass):
    n = prob_mass.shape[0]
    result = np.zeros(n)
    for i, e in enumerate(prob_mass):
        result += np.concatenate([np.zeros(i), e*np.ones(n-i)])
    #result -= result.min()
    return result
def pmf2cdf_2d(prob_mass):
    n, m = prob_mass.shape
    result = np.zeros(shape=(n, m))
    for i in range(0, n):
        for j in range(0,m):
            row = np.zeros(shape=(i, m)) if i != 0 else None
            col = np.zeros(shape=(n-i, j)) if j !=0 else None

            add_mass = prob_mass[i,j] * np.ones(shape=(n-i, m-j))

            if col is not None:
                add_mass = np.concatenate([col, add_mass], axis=1)
                
            if row is not None:
                add_mass = np.concatenate([row, add_mass], axis=0)
            
            result += add_mass
    result -= result.min()
    return result

def pmf_cond(i, pro_mass, axis =0):# axis=0: x, axis=1: y
    n,m = pro_mass.shape
    if axis==0:
        if not(i>=0 and i<n):
            raise IndexError("Exceeded data dimension.")
        pdf = pro_mass[i]
    elif axis ==1:
        if not(i>=0 and i<m):
            raise IndexError("Exceeded data dimension.")
        pdf = pro_mass[:, i]
    
    return pdf /pdf.sum()

def cdf_cond(i, pro_mass, axis=0):# axis=0: x, axis=1: y
    return pmf2cdf(pmf_cond(i, pro_mass, axis))


# Chebyshev approximation
def get_cheb_approx_pdf(pos_mass, x, region=(None, None)):
    if pos_mass.shape != x.shape:
        raise ValueError("Dimensions are not same each others.")
    
def get_cheb_approx_cdf_pdf(pos_mass, x, region=(None, None)):
    if pos_mass.shape != x.shape:
        raise ValueError("Dimensions are not same each others.")    
    
def get_cheb_approx_cdf_pmf(pos_mass, dx, region=(None, None)):
    return 0


# Inverse Transform Sampling

# Direct
def int_sampling(
        uni_sam, 
        pmf= None,
        cdf = None,
        domain = None,
        interpolate = False # 1 dim interpolate between sample points
        ):
    if pmf.min() <0 or pmf.max() >1:
        raise ValueError("Invaild value in probability mass vector. Exceeding range [0,1]")
    if math.fabs(pmf.sum() -1) > EPS:
        raise ValueError("Invaild probability mass vector. The sum is not 1.")

    if uni_sam.min() <0 or uni_sam.max() >1:
        raise ValueError("Invaild value in sample data. Exceeding range [0,1]")
    if pmf is None and cdf is None:
        raise ValueError("")
    
    if isinstance(domain, Iterable) and len(domain) ==2:
        domain = np.linspace(domain[0], domain[1], endpoint=True)

    if cdf is None:
        cdf = pmf2cdf(pmf)
    
    if domain is not None and domain.size != cdf.size:
        raise ValueError("Dimension error, domain and cdf have different dimension.")
        
    n = cdf.size

    if interpolate:
        samples = []
        for uni_i in uni_sam:
            if uni_i <= cdf.min():
                samples.append(0)
                continue
            if uni_i >= cdf.max():
                samples.append(n-1)
                continue
            
            min_sol = np.where(cdf <=uni_i)[0]
            max_sol = np.where(cdf >=uni_i)[0]


            min_index = min_sol.max() if min_sol.size != 0 else None
            max_index = max_sol.min() if max_sol.size != 0 else None

            if min_index == max_index:
                if min_index is None:
                    pass # impossible
                else:
                    samples.append(min_index)
            elif min_index is None or max_index is None:
                if min_index is None:
                    samples.append(max_index)
                else:
                    samples.append(min_index)
            else:
                dy = cdf[max_index] - cdf[min_index]
                dy_i = uni_i - cdf[min_index]
                dx_i = dy_i / dy
                samples.append(min_index + dx_i)
        result = np.array(samples)
    else:
        results = []
        uni_sam.sort()

        for i, c_i in enumerate(cdf):
            uni_index = np.where(uni_sam <= c_i)[0]
            if uni_index.size == 0:
                continue
            results.append(i*np.ones(uni_index.size))
            uni_sam = uni_sam[uni_index.size:]

        if domain is not None:
            for i in range(0, len(results)):
                results[i] *= domain[i]

        result = np.concatenate(results)
        
    return result

# Continuous version
## Root finding
def int_sampling_root(
        uni_sam,
        cdf:Callable,
        domain
    ):
    a, b = domain
    if a > b:
        a, b = b, a

    cdf = lambda x: cdf()
    result = []
    for uni_s in uni_sam:
        result.append(bisect(lambda x: cdf(x) - uni_s, *domain))
    return np.array(result)

# SVD Chebyshev approximation
#
SingleFunction = Callable[[Number], Number]
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
        fx_v = np.array([list(map(f, x)) for f in self._fx_list])
        fy_v = np.array([list(map(f, y)) for f in self._fy_list])
        return (weight * fx_v * fy_v).sum(axis=0).reshape(shape)
        
    def __check_types(self, fx_list:NDArray, fy_list:NDArray, weights:NDArray):
        for fx in fx_list:
            if isinstance(fx, Callable):
                continue
            raise TypeError("fx list must consist of callable objects.")
        for fy in fy_list:
            if isinstance(fy, Callable):
                continue
            raise TypeError("fy list must consist of callable objects.")
        for w in weights:
            if isinstance(w, Number):
                continue
            raise TypeError("weight must consist of number type.")
    @classmethod
    def from_pivots(cls, xy_pivots, f:Callable, domain, cheby_deg):
        return cls( domain=domain, cheby_deg = cheby_deg)
    def rank_up(self, 
                        fx_list:Union[SingleFunction, Iterable[SingleFunction]], 
                        fy_list:Union[SingleFunction, Iterable[SingleFunction]], 
                        weights:Union[Number, Iterable[Number]]):
        fx_list = list(fx_list) if isinstance(fx_list, Iterable) else [fx_list]
        fy_list = list(fy_list) if isinstance(fy_list, Iterable) else [fy_list]
        weights = list(weights) if isinstance(weights) else [weights]
        if (fx_list.shape) != 1:
            fx_list = fx_list.reshape(-1) 
        if (fy_list.shape) != 1:
            fy_list = fy_list.reshape(-1) 
        if (weights.shape) != 1:
            weights = weights.reshape(-1) 
        l_a = len(fx_list) 
        l_b = len(fy_list)
        l_c = len(weights)

        if l_a != l_b or l_b != l_c or l_a != l_c:
            raise ValueError("Dimensions are not same, {l_a}, {l_b}, {l_c}.")
        
        self.__check_types(fx_list, fy_list, weights)

        self._fx_list += np.concatenate([self._fx_list, fx_list])
        self._fy_list += np.concatenate([self._fy_list, fy_list])
        self._weights += np.concatenate([self._weights, weights])
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
    

def function_vec(f_arr, x, multi=False): # Array of functions apply to x value
    if multi:
        return np.fromiter((fi(*x) for fi in f_arr), dtype=np.dtype(type(x)))
    else:
        return np.fromiter((fi(x) for fi in f_arr), dtype=np.dtype(type(x)))
def cal_f_xy_approx(x, y, f_x, f_y, coefs):
    f_x_v = function_vec(f_x, x)
    f_y_v = function_vec(f_y, y)
    return coefs.dot(f_x_v.dot(f_y_v))
f_xy_approx_vec = np.vectorize(cal_f_xy_approx, excluded=["f_x", "f_y", "coefs"])

def u(x, f, domain, deg): # return u(y)=f(x_i, y) Chebyshev approximation
    a, b, c, d = domain
    
    y_list = np.linspace(b, d, deg+1, endpoint=True)
    x_list = np.full_like(y_list, x)
    z_list = f(np.hstack(x_list, y_list))
    u_i = Chebyshev.fit(y_list, z_list, deg=deg, domain=domain)
    return u_i
def v(y, f, domain, deg):# return v(x)= f(x, y_i) Chebyshev approximation
    a, b, c, d = domain
    
    x_list = np.linspace(a, c, deg+1, endpoint=True)
    y_list = np.full_like(x_list, y)
    z_list = f(np.hstack(x_list, y_list))
    v_i = Chebyshev.fit(x_list, z_list, deg=deg, domain=domain)
    return v_i
def u_i(xy_list, index, f, domain, deg):
    x, y = xy_list[index]
    return u(x, f, domain, deg)
def v_i(xy_list, index, f, domain, deg):
    x, y = xy_list[index]
    return v(y, f, domain, deg)
def ge_approx(xy_points, f, domain, chebyshev_deg=100):
    z = f(*xy_points.transpose())
    max_i = np.argmax(z)
    u_func = [u_i(xy_points, max_i, f, domain, chebyshev_deg)]
    v_func = [v_i(xy_points, max_i, f, domain, chebyshev_deg)]
    d_coef = [1/(f(*xy_points[max_i]))]
    for point in xy_points:
        pass


