from typing import Callable
from collections.abc import Iterable
import numpy as np


# LED irradiation models 
LED_VERTICAL_AXIS = np.array([0,0,1])
def cosine(target:np.array):
    r = target
    r2 = np.dot(r,r)
    r_size = np.sqrt(r2)
    return np.dot(r/r_size, LED_VERTICAL_AXIS)
def theta(target:np.array):
    return np.arccos(cosine(target))
def lambertian(target:np.array,s:np.double):
        return np.power(cosine(target), s)
def gaussian(target:np.array, h:np.double):
        theta = theta(target)
        return np.exp(-np.power(theta/h,2))

class LED: #Irradiation model of LED 
    def __init__(self, I:float, location:np.ndarray, direction:np.ndarray, model_name:str, model:Callable, **model_params):
        self.I = I
        self._location = location
        self._direction = direction
        self._rotation_rel = self._inverse_rotation(self._direction)

        self._model_name = model_name
        self._model = model
        self._model_params = model_params
    def intensity(self, target: np.ndarray):
        #Convert to relative vectors for model
        r_rel = self._rotation_rel@(target - self._location)
        r2 = r_rel@r_rel
        if len(self._model_params) == 0:
            intense = self.I/r2 * self._model(r_rel)
        else:
            intense = self.I/r2 * self._model(r_rel, self._model_params)
        return intense
    def _inverse_rotation(self, direction):
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
        return self._location
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
        self._location = new_location
    @property
    def direction(self):
        return self._direction
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
        self._direction = new_direction/d
        self._rotation_rel = self._inverse_rotation(self._direction)
    @property
    def model_info(self):
        info ={
            "name": self._model_name,
            "parameters": self._model_params
        }
        return info
    @property
    def model_parameters(self):
        return self._model_params
    @model_parameters.setter
    def model_parameters(self, **new_params):
        if len(self._model_params) != len(new_params):
            raise ValueError("Parameters number is not same")
        new_keys = new_params.keys()
        old_keys = self._model_params.keys()
        for key in new_keys:
            if key not in old_keys:
                raise KeyError("Invaild key exist in new parameters.")
            self._model_params[key] = new_params[key]

    @classmethod
    def Lambertian(cls, I, s, location, direction=LED_VERTICAL_AXIS):
        model_type = "Lambertian"
        return cls(I, location, direction, model_type, lambertian, s=s)
    @classmethod
    def Gaussian(cls, I, h, location, direction=LED_VERTICAL_AXIS):
        model_type = "Gaussian"
        return cls(I, location, direction, model_type, gaussian, h=h)

class LEDarray:
    def __init__(self, array, location=np.array([0,0,0]), direction = np.array([0,0,1])):
        self.array = array
        self.N = array.shape[1]
        self.M = array.shape[0]
        
        self.location = location
        self.direction = direction

    @classmethod
    def fromAxisArray(cls, I0, xarray, yarray=np.array([0.0])):
        array = np.array([[[i, j] for i in xarray] for j in yarray])
        return cls(s, I0, array)
    def ledcoordinate(self, i, j):
        return self.array[j,i]

    