import math
from typing import Callable
from collections.abc import Iterable

import numpy as np


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
# Root finding
from scipy.optimize import bisect
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
    
def int_sampling_2d(uni_sam, cdf, order=0):#
    return 0
#    Chebyshev
def int_cheby_sampling():
    return 0
def int_cheby_sampling_2d():
    return 0

# Continuous version


# SVD Chebyshev approximation

class RankApprox2dim:
    def __init__(self, fx_list, fy_list, weights):
        self.fx_list = fx_list
        self.fy_list = fy_list
        self.weights = weights
        pass
    @classmethod
    def from_pivots(cls, xy_pivots, f:Callable, domain, cheby_deg):
        return cls()
    
    def __check_types(self, fx_list, fy_list, weights):
        for fx in fx_list:
            if isinstance(fx, Callable):
                continue
            raise TypeError("fx list must consist of callable objects.")
        for fy in fy_list:
            if isinstance(fy, Callable):
                continue
            raise TypeError("fy list must consist of callable objects.")
        for w in weights:
            if isinstance(w, Numeric)
        pass
    def add_elements(self, fx_list, fy_list, weights):
        fx_list = list(fx_list) if isinstance(fx_list, Iterable) else [fx_list]
        fy_list = list(fy_list) if isinstance(fy_list, Iterable) else [fy_list]
        weights = list(weights) if isinstance(weights) else [weights]

        l_a = len(fx_list) 
        l_b = len(fy_list)
        l_c = len(weights)

        if l_a != l_b or l_b != l_c or l_a != l_c:
            raise ValueError("Dimensions are not same, {l_a}, {l_b}, {l_c}.")
        
        self.__check_types(fx_list, fy_list, weights)

        self.fx_list += fx_list
        self.fy_list += fy_list
        self.weights += weights

    def fx(self, x):
        self.fx_list
        return result
    def fy(self, y):
        self.fy_list
        return result
    def __call__(self, x, y):
        f_x_v = self.f_x

        return self.coefs.dot( self.f_x.dot())



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


