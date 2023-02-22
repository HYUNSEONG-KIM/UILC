from typing import Tuple, Literal, Callable

import numpy as np
from numpy.lib.stride_tricks import as_strided

# Matrix
def toeplitz(c, r):
    # Imgrated source code from scipy `toeplitz` function
    # BSD-licensed
    c = np.asarray(c).ravel()
    if r is None:
        r = c.conjugate()
    else:
        r = np.asarray(r).ravel()
    # Form a 1-D array containing a reversed c followed by r[1:] that could be
    # strided to give us toeplitz matrix.
    vals = np.concatenate((c[::-1], r[1:]))
    out_shp = len(c), len(r)
    n = vals.strides[0]
    return as_strided(vals[len(c)-1:], shape=out_shp, strides=(-n, n)).copy()

# Utils
def restore_data(data, dim_k, dim_crop):
    l, k = dim_k
    cr, cc = dim_crop
    er = l - cr
    ec = k -cc
    return data[er: -er, ec: -ec]

def get_dim_ext(dim_k, dim_crop):
    l, k = dim_k
    cr, cc = dim_crop

    if cr * cc *l *k  == 0:
        raise ValueError(f"None of dimension can be zero, l:{l}, k:{k}, cc:{cc}, cr:{cr}") 
    if cr <0 or cc<0 or l<0 or k<0:
        raise ValueError(f"None of dimension can be negative, l:{l}, k:{k}, cc:{cc}, cr:{cr}") 
    if cr > l or cc >k:
        raise ValueError(f"Cropping dimension must be smaller than kernel dimension.")
    return l-cr, k-cc

def get_dim_convole(dim_i, dim_k, dim_c):
    n, m = dim_i
    l, k = dim_k

    if l> n or k>m:
        n ,l = l, n
        m, k = k, m

    cr, cc = dim_c
    
    return n+l - 2*cr+1, m+k - 2*cc +1

# Internal routines
def _expand_matrix(data, dim_ext, edge_param):
    er, ec = dim_ext

    if edge_param[0] == "custom" and len(edge_param) >=2:
        C1, C2, R1, R2  = edge_param[1](data, (er, ec), edge_param[2:])
    else: C1, C2, R1, R2 = _edge_matrix(data, (er, ec), edge_param)

    # Column expand
    if isinstance(C1, np.ndarray):
        data = np.concatenate([C1, data, C2], axis=1)
    # Row expand
    if isinstance(R1, np.ndarray):
        data = np.concatenate([R1, data, R2], axis=0)

    return data

def _edge_matrix(data, dim_ext, edge_param):
    n, m = data.shape
    er, ec = dim_ext
    # Column Matrix: (n, ec)
    # Row Matrix: (er, m+2ec) 
    dim_column = (n, ec)
    dim_row = (er, int(m+2*ec))

    e_method_name, e_method_param = edge_param

    if e_method_name == "extend":
        
        a_col_1 = np.expand_dims(data[0:, 0], axis=0)
        a_col_2 = np.expand_dims(data[0:,-1], axis=0)
        C1 = np.repeat(a_col_1, ec, axis=0).transpose()
        C2 = np.repeat(a_col_2, ec, axis=0).transpose()

        a_row_1 = np.expand_dims(np.concatenate([C1[0],data[0], C2[0]]), axis=0)
        a_row_2 = np.expand_dims(np.concatenate([C1[-1], data[-1], C2[-1]]), axis=0)
        R1 =np.repeat(a_row_1, er, axis=0)
        R2 =np.repeat(a_row_2, er, axis=0) 

    elif e_method_name ==  "wrap":
        if er >0 and ec >0:
            Cor1 = data[-er:,   -ec:  ] 
            Cor2 = data[-er:,      :ec] 
            Cor3 = data[   :er, -ec:  ]
            Cor4 = data[   :er,    :ec]

            Row1 = data[-er:,    :]
            Row2 = data[   :er,  :]

            R1 = np.concatenate([Cor1, Row1, Cor2], axis=1)
            R2 = np.concatenate([Cor3, Row2, Cor4], axis=1)

            C1 = data[:, -ec:  ] 
            C2 = data[:,    :ec]
        elif er>0: # ec ==0
            R1 = data[-er:,    :]
            R2 = data[   :er,  :]
            C1 = C2 = None
        else:
            R1 = R2 = None
            C1 = data[:, -ec:  ] 
            C2 = data[:,    :ec]

    elif e_method_name == "reflect":
        if er * ec >0:
            C1 = np.flip(data[:,    :ec], axis=1)
            C2 = np.flip(data[:, -ec:  ], axis=1)

            Cor1 = np.flip(C1[   :er, :], axis=0)
            Cor2 = np.flip(C2[   :er, :], axis=0)
            Cor3 = np.flip(C1[-er:  , :], axis=0)
            Cor4 = np.flip(C2[-er:  , :], axis=0)

            Row1 = np.flip(data[:er, :], axis=0)
            Row2 = np.flip(data[-er:, :], axis=0)

            R1 = np.concatenate([Cor1 , Row1, Cor3], axis =1) 
            R2 = np.concatenate([Cor2 , Row2, Cor4], axis =1)
        elif er > 0:
            R1 = np.flip(data[:er, :], axis=0)
            R2 = np.flip(data[-er:, :], axis=0)
            C1 = C2 = None
        else:
            R1 = R2 = None
            C1 = np.flip(data[:,    :ec], axis=1)
            C2 = np.flip(data[:, -ec:  ], axis=1)
    elif e_method_name == "mirror":
        if er*ec >0:
            C1 = np.flip(data[:,    1:ec+1], axis=1)
            C2 = np.flip(data[:, -(ec+1): -1  ], axis=1)

            Cor1 = np.flip(C1[   :er, :], axis=0)
            Cor2 = np.flip(C2[   :er, :], axis=0)
            Cor3 = np.flip(C1[-er:  , :], axis=0)
            Cor4 = np.flip(C2[-er:  , :], axis=0)

            Row1 = np.flip(data[1:er+1, :], axis=0)
            Row2 = np.flip(data[-(er+1):-1, :], axis=0)

            R1 = np.concatenate([Cor1 , Row1, Cor3], axis =1) 
            R2 = np.concatenate([Cor2 , Row2, Cor4], axis =1)
        elif er >0:
            R1 = np.flip(data[1:er+1, :], axis=0)
            R2 = np.flip(data[-(er+1):-1, :], axis=0)
            C1 = C2 = None
        else:
            R1 = R2 = None
            C1 = np.flip(data[:,    1:ec+1], axis=1)
            C2 = np.flip(data[:, -(ec+1): -1  ], axis=1)
    elif e_method_name == "constant":
        c = e_method_param
        C1 = C2 = c*np.ones(shape=dim_column) if ec>0 else None
        R1 = R2 = c*np.ones(shape=dim_row) if er>0 else None
    
    else:
        C1 = C2 = R1 = R2 = None

    return C1, C2, R1, R2

def convolve2d(  
    data:np.ndarray, 
    filter:np.ndarray,
    crop:Tuple[int, int]=(1, 1),
    edge_mode:Literal["extend", "wrap", "mirror","constant"]="constant",
    edge_params = [0]
    ):
    l, k = filter.shape
    er, ec = get_dim_ext((l, k), crop)

    data_ext = _expand_matrix(data, (er, ec), [edge_mode]+edge_params)
    n,m = data_ext.shape

    result = np.zeros(shape = (n-l+1, m-k+1))
    for i, h_row in enumerate(filter):
        data_i = data_ext[i:n-l+i+1]
        result += np.stack([np.convolve(r, h_row, mode="vaild") for r in data_i])

    return result

#

def _vec2sub_toeplitz(vec, m): #m: column dimension of matrix
    k = len(vec)
    n = m-k +1
    rows = []
    for i in range(0, n):
        row_zeros =  np.zeros(shape=(m-k,))
        rows.append(np.insert(row_zeros, i, vec))
    return np.vstack(rows)


def convolve2toeplitz(
    data:np.ndarray, 
    filter:np.ndarray, 
    crop:Tuple[int, int]=(1, 1),
    edge_mode:Literal["extend", "wrap", "mirror","constant"]="constant",
    edge_params = [0], 
    preserve_filter=False) -> Tuple[np.ndarray, np.ndarray, Callable]:
    
    l, k = filter.shape
    er, ec = get_dim_ext((l, k), crop)
    data_ext = _expand_matrix(data, (er, ec), [edge_mode]+edge_params)

    n,m = data_ext.shape

    H_list = []
    for i in range(0, l):
        H_list.append(_vec2sub_toeplitz(filter[i], m))
    
    topelitz_dim = H_list[0].shape

    rows =[]
    for i in range(0, n-l+1):
        i_f = i
        i_b = n- l -i_f

        row = i_f*[np.zeros(shape=topelitz_dim)] + H_list + i_b*[np.zeros(shape=topelitz_dim)]
        rows.append(row)
    return np.block(rows), data_ext.flatten()

#-------------------------------------------------------------------------------------------------------
# Special Case
#---------------------------------------------------------------

def get_matrix_system(filter, dim_d):
    n,m = dim_d
    l, k= filter.shape

    if k != 2*m-1 or l != 2*n-1:
        raise ValueError("Invaild dimension: l, k must be 2n-1, 2m-1")
    
    rows = []
    for i in range(0, n):
        row_i = n-1-i
        row_f = 2*n-1-i

        for j in range(0, m):
            column_i = m-1 -j
            column_f = 2*m-1 -j

            #print(row_i, row_f)
            #print(column_i, column_f)

            t = filter[row_i : row_f, column_i:column_f]

            #print(t.shape)
            rows.append(t.flatten())

    #mat = np.vstack(rows)
    return np.vstack(rows)

