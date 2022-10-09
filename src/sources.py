from typing import Callable
from collections.abc import Iterable
import numpy as np

class LED: #Irradiation model of LED 
    def __init__(self, I:float, location:np.ndarray, direction:np.ndarray, model_name:str, model:Callable, **model_params):
        self.I = I
        self.__location = location
        self.__direction = direction
        self.__rotation_rel = self.__inverse_rotation(self.__direction)

        self.__model_name = model_name
        self.__model = model
        self.__model_params = model_params
    def intensity(self, target: np.ndarray):
        #Convert to relative vectors for model
        r_rel = self.__rotation_rel@(target - self.__location)
        r2 = r_rel@r_rel
        return self.I/r2 * self.__model(r_rel, self.__model_params)
    def __inverse_rotation(self, direction):
        x, y, z = direction
        theta_xy = np.arctan(y/x)
        theta_zx = np.arccos(z)

        c_z = np.cos(-theta_xy)
        s_z = np.sin(-theta_xy)
        c_y = np.cos(-theta_zx)
        s_y = np.sin(-theta_zx)

        R_y_inverse = np.array(
            [[c_y, 0, s_y],
             [ 0 , 1,  0 ],
             [-s_y, 0, c_y]])
        R_z_inverse = np.array(
            [[c_z, -s_z, 0],
             [s_z, c_z,  0],
             [0,    0,   1]])
        return R_y_inverse@R_z_inverse
    @property
    def location(self):
        return self.__location
    @location.setter
    def location(self, new_location):
        if isinstance(new_location) != np.ndarray:
            if isinstance(new_location, Iterable):
                try:
                    new_location = np.array(new_location)
                except:
                    raise TypeError("Cannot convert to numpy array")
        if len(new_location) != 3:
            raise ValueError(f"Length must be 3, current: {len(new_location)}")
        self.__location = new_location
    @property
    def direction(self):
        return self.__direction
    @direction.setter
    def direction(self, new_direction):
        if isinstance(new_direction) != np.ndarray:
            if isinstance(new_direction, Iterable):
                try:
                    new_direction = np.array(new_direction)
                except:
                    raise TypeError("Cannot convert to numpy array")
        if len(new_direction) != 3:
            raise ValueError(f"Length must be 3, current: {len(new_direction)}")
        d = np.sqrt(np.dot(new_direction, new_direction))
        self.__direction = new_direction/d
        self.__rotation_rel = self.__inverse_rotation(self.__direction)
    @property
    def model_info(self):
        info ={
            "name": self.__model_name,
            "parameters": self.__model_params
        }
        return info
    @property
    def model_parameters(self):
        return self.__model_params
    @model_parameters.setter
    def model_parameters(self, **new_params):
        if len(self.__model_params) != len(new_params):
            raise ValueError("Parameters number is not same")
        new_keys = new_params.keys()
        old_keys = self.__model_params.keys()
        for key in new_keys:
            if key not in old_keys:
                raise KeyError("Invaild key exist in new parameters.")
            self.__model_params[key] = new_params[key]
    

    @classmethod
    def Lambertian(cls, I, s, location, direction=np.array([0,0,1])):
        model_type = "Lambertian"
        def lambertian(s:np.double, target:np.array):
            r = target
            r2 = np.dot(r,r)
            r_size = np.sqrt(r2)
            cos = np.dot(r/r_size, np.array([0,0,1]))
            return np.power(cos, s)
        return cls(I, location, direction, model_type, lambertian, s)
    #@classmethod
    #def Gaussian(cls, I, location, direction=np.array([0, 0 ,1])):
    #    model_type = "Gaussian"
    #    def gaussian()
    #    return cls(I, location, direction, model_type, )

class LEDarray:
    def __init__(self, s, I0, array):
        self.array = array
        self.N = array.shape[1]
        self.M = array.shape[0]
        self.s = s
        self.I0 = I0
    @classmethod
    def fromAxisArray(cls, s, I0, xarray, yarray=np.array([0.0])):
        array = np.array([[[i, j] for i in xarray] for j in yarray])
        return cls(s, I0, array)
    def ledcoordinate(self, i, j):
        return self.array[j,i]

    