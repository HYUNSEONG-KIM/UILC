import numpy as np
import scipy as scp
from admath import ASM, hyperg_2F1


class Hyeon:
    def __init__(self, led, width, height=200E-3):
        self.case = case
        self.led = led
        self.nled = nled
        
        self.width =width
        self.height =height
        self.ledarray = self.calarray()

    def beta(self, x):
        k = self.width/self.height
        d = x/self.height
        p = -(self.led.m/2 +1)

        return 1/2 *(pow((1+d**2)/(1+ (k/2 -d)**2), p ) + pow((1+ d**2)/(1+(k/2 +d)**2), p))
    
    def get_ds(self):
        
        sol = scp.optimizae.root_scalar(
            lambda x: self.beta(x)-1, 
            bracket = [0, self.width*0.7], method="brentq")
        return sol.root 
    
    def calarray(self):
        

    def intensity(self, target):
        print()
    def grad_intensity(self, target, axis=0, ):