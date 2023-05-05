import math
from functools import partial

from typing import Callable, Tuple, Union
from numpy.typing import NDArray, ArrayLike
from numbers import Number
from collections.abc import Iterable

import numpy as np

from numpy.polynomial.chebyshev import Chebyshev
from scipy.optimize import bisect
from convolution import convolve2d

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
    if pmf is not None:
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

    if cdf is None and pmf is not None:
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
    a, b = domain.min(), domain.max()
    width = b-a
    result = result/width -a
    return result

# Continuous version
## Root finding
def int_sampling_root(
        uni_sam,
        cdf:Callable,
        domain= [-1, 1]
    ):
    a, b = domain
    if a > b:
        a, b = b, a

    #cdf = lambda x: cdf(x)
    result = []
    for uni_s in uni_sam:
        try:
            sol = bisect(lambda x: cdf(x) - uni_s, *domain)
        except ValueError:
            result.append(a) if uni_s <0.5 else result.append(b)
        result.append(bisect(lambda x: cdf(x) - uni_s, *domain))
    return np.array(result)
