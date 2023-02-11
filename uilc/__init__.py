import math
from typing import Tuple, Literal, Union
from collections.abc import Iterable

import numpy as np
from scipy.optimize import root_scalar, minimize_scalar

from uilc.utils.misc import float_eps, csym_index

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
        
        xarr = dx*(csym_index(N))
        yarr = dy*(csym_index(M))
        
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


__all__ =[
    "PositionArray",
    "utils",
    "crit",
    "radiation"
]