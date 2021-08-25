import numpy as np

#---------------------------------------LED
class Lamber:
    def __init__(self, intensity, m, location):
        self.intense  = intensity
        self.m = m
        if location.size == 3:
            self.location = location #3 dimenstional array 

    def recalltype(self):
        return "Lamber"

    def setlocation(self, location):
        if location.nsize == 3:
            self.location = location #3 dimenstional array 
    
    def intensity(self, theta):
        if theta < 0 or theta > np.pi/2:
            return "range error"

        return self.intense * np.power(np.cos(theta),self.m)

    def target_intensity(self, target):
        if target.size != 3:
            return "err"
        direction = target - self.location
        r2 = np.dot(direction,direction)
        theta = np.arccos(direction[2]/np.sqrt(r2))
        return self.intensity(theta) / r2

class Poly:
    def __init(self, intensity, a_n , location):
        self.intense = intensity
        if a_n.ndim !=1:
            return "err"
        self.n = a_n.size
        self.a_n = a_n
        if location.size == 3:
            self.location = location #3 dimenstional array
    
    def recalltype(self):
        return "Poly"
    
    def setlocation(self, location):
        if location.size == 3:
            self.location = location #3 dimenstional array
    
    def intensity(self,theta):
        if theta < 0 or theta > np.pi/2:
            return "range error"
        
        cof = self.a_n[self.n]
        for i in range(self.n,0,1):
            cof = cof * theta + self.a_n[i]
        return self.intense * cof
    
    def target_intensity(self, target):
        if target.size != 3:
            return "err"
        direction = target - self.location
        r2 = np.dot(direction,direction)
        theta = np.arccos(direction[2]/np.sqrt(r2))
        return self.intensity(theta) / r2
 