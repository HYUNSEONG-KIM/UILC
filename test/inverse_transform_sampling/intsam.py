import math

import numpy as np


EPS = 10E6*np.finfo(float).eps

# Discrete version
def pmf2cdf(prob_mass):
    n = prob_mass.shape[0]
    result = np.zeros(n)
    for i, e in enumerate(prob_mass):
        result += np.concatenate([np.zeros(i), e*np.ones(n-i)])
    return result/result.max()
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
    return result/result.max()

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

#    Direct
def int_sampling(uni_sam, pos_mass):
    if pos_mass.min() <0 or pos_mass.max() >1:
        raise ValueError("Invaild value in probability mass vector. Exceeding range [0,1]")
    if math.fabs(pos_mass.sum() -1) > EPS:
        raise ValueError("Invaild probability mass vector. The sum is not 1.")

    if uni_sam.min() <0 or uni_sam.max() >1:
        raise ValueError("Invaild value in sample data. Exceeding range [0,1]")
    
    cdf = pmf2cdf(pos_mass)
    cdf /= cdf.max()
    n = len(cdf)

    samples = []

    for uni_i in uni_sam:
        if uni_i <= cdf.min():
            samples.append(0)
        if uni_i >= cdf.max():
            samples.append(n-1)
        
        min_sol = np.where(cdf <=uni_i)[0]
        max_sol = np.where(cdf >=uni_i)[0]

        
        min_index = min_sol.max() if len(min_sol) >0 else None
        max_index = max_sol.min() if len(max_sol) >0 else None

        
        
        if min_index == max_index:
            if min_index is None:
                pass
            else:
                samples.append(min_index)
        elif min_index is None or max_index is None:
            if min_index is None:
                samples.append(max_index)
            else:
                samples.append(min_index)
        else:
            eps_min, eps_max = math.fabs(cdf[min_index]-uni_i), math.fabs(cdf[max_index]-uni_i)

            if eps_min < eps_max:
                samples.append(min_index)
            else:
                samples.append(max_index)
        
    return np.array(samples)
def int_sampling_2d():
    return 0
#    Chebyshev
def int_cheby_sampling():
    return 0
def int_cheby_sampling_2d():
    return 0

# Continuous version