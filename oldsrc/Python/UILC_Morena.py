from UILC_general import *


class Morena:
    self.L = "Linear"
    self.R = "Rectangle"

    def __init__(self, led, nled, case="Linear"):
        self.case = case
        self.led = led
        self.nled = nled
        self.dm = self.dm()
        self.array = self.array()
    
    def dm(self, deep):
        if self.case == self.L:
            #get linear dm
            return dm
        elif self.case == self.R:
            #get Rectangle dm\
            return dm

    def modify_lednum(self, n):
        self.nled = n
        self.dm = self.dm()
        self.array = self.array()
    
    def get_area(self):
        if self.case == self.L:
            #get linear dm
            return self.dm * (self.nled -1)
        elif self.case == self.R:

    def array(self,height):
    
    def intensity(self, target):
        #calculate the intensity of the array 
        I=0
        for i in range(0,len(self.array)):
            I = I + self.led.intensity(target)

        return I 
    
    def grad_intensity(self, target, axis=0):
        if axis ==0: #x axis
        elif axis == 1: #y axis
        elif axis ==2: #xy gradient
        else:
            return "err"
    
    
        
