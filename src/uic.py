import math
import numpy as np
from scipy import optimize as op

# This LED class represent imperpect Lambertian model with inverse-square law
# For hemosphere lens and far field irradiation case, it fits well such situations.
# the theta, phi represent the orientation of LED, normally it is considered as positive z axis direction.
# theta: angle between xy projection vector of LED orientation and x-axis. Range= [0-2pi)
# phi: angle between z-axis and LED orient vector. Range = [0,pi/2)
class LED:
    def __init__(self, I, m, location, theta, phi):
        if isinstance(I,(int,float)) and isinstance(m,(int,float)) and isinstance(phi,(int,float)):
            if I >0 and m >0 and abs(phi) < math.pi/2:
                pass
            else:
                raise ValueError("The argument I, m must be positive value and phi is an angle lying on [-pi/2,pi/2]\n I={}, m ={}, phi={}".format(I,m,phi))
        else:
            raise ValueError("The argument I,m,phi must be numerical type int, float:\n I={}, m ={}, phi={}".format(type(I),type(m),type(phi)))

        if type(location) is not np.ndarray:
            raise ValueError("The led location must be numpy array ndarray with length 3.\n location variable type:{}".format(type(location)))
        elif location.size !=3:
            raise ValueError("The vector size must be 3.\n ndarray size:{}".format(location.size))
        elif not isinstance(location[0],(int,float,np.integer)):
            raise ValueError("The vector type must be numerical type float or int.\n ndarray type:{}, element type:{}".format(location.dtype,type(location[0])))
        
        self.I=I
        self.m =m 
        self.location = location 
        self.phi = phi
        self.theta =theta
        self.ori = np.array([math.sin(self.phi)*math.cos(self.theta), math.sin(self.phi)*math.sin(self.theta), math.cos(self.phi)])
    
    def intensity(self, target):
        # Value checker
        if type(target) is not np.ndarray:
            raise ValueError("The target argument must be numpy array type.")
        elif target.size !=3:
            raise ValueError("The vector size must be 3.\n ndarray size:{}".format(target.size))
        elif not isinstance(target[0],(int,float,np.integer)):
            raise ValueError("The vector type must be numerical type float or int.\n ndarray type:{}, element type:{}".format(target.dtype,type(target[0])))

        d= target - self.location
        r2 = np.dot(d,d)
        dsize=math.sqrt(r2)
        cos = np.dot(d,self.ori)/dsize

        return self.I/r2 * math.pow(cos,self.m)

    def set_orientation_angle(self, theta, phi):
        # Value checker

        self.theta = theta
        self.phi = phi
        self.ori = np.array([math.sin(self.phi)*math.cos(self.theta), math.sin(self.phi)*math.sin(self.theta), math.cos(self.phi)])

    def set_orientation_coors(self, orient):
        # Value checker
        if type(orient) is not np.ndarray:
            raise ValueError("The led orient must be numpy array ndarray with length 3.\n location variable type:{}".format(type(orient)))
        elif orient.size !=3:
            raise ValueError("The vector size must be 3.\n ndarray size:{}".format(orient.size))
        elif not isinstance(orient[0],(int,float,np.integer)):
            raise ValueError("The vector type must be numerical type float or int.\n ndarray type:{}, element type:{}".format(orient.dtype,type(orient[0])))  
        
        self.orient = orient/math.sqrt(np.dot(orient,orient))
        self.phi = math.acos(self.ori[2])
        self.theta = math.acos(self.ori[0]/math.sqrt(1-self.ori[2]**2))

class LEDmatrix:
    def __init__(self, m, I0,  xarray, yarray=np.array([[0]])):
        self.N = xarray.shape[1] 
        self.M = yarray.shape[1]
        self.lednum = self.N*self.M 
        self.m = m
        self.matrix, self.dim = self.initLEDArray(m,I0, xarray, yarray) # 2dim ndarray
        self.area = self.dim[0] * self.dim[1] 
    
    def initLEDArray(self, m, I, xarr,yarr):
        ledmatrix = []
        for j in yarr[0]:
            ledarray=[]
            for i in xarr[0]:
                ledarray.append(LED(I, m, np.array([i,j,0]), 0, 0))
            ledmatrix.append(ledarray)
        vec = ledmatrix[yarr.shape[1]-1][xarr.shape[1]-1].location - ledmatrix[0][0].location

        return ledmatrix, (vec[0],vec[1])

    def LED_location(self, i, j):
        return self.matrix[i][j].location

    def LED_graphic_view(self):
        pass
        
    def intensity(self, target, position=np.array([0,0,0]), orient = np.array([0,0,1])):
        r = target - position
        orient = orient/math.sqrt(np.dot(orient,orient))
        phi = math.acos(orient[2])
        if math.isclose(phi,0,rel_tol=np.finfo(float).eps):
            theta =0
        else:
            theta = math.acos(orient[0] / math.sqrt(orient[0]**2 + orient[1]**2))
            theta = theta if orient[1] >0 else 2*math.pi - theta

        sp = math.sin(phi)
        cp = math.cos(phi)
        st = math.sin(theta)
        ct = math.cos(theta)

        Ryp = np.array([[cp , 0,sp],
                        [0  , 1, 0],
                        [-sp, 0,cp]])

        Rzt = np.array([[ct,-st,0],
                        [st, ct,0],
                        [0 ,  0,1]])

        x1 = Rzt.dot(Ryp.dot(np.array([1,0,0])))
        y1 = Rzt.dot(Ryp.dot(np.array([0,1,0])))
        z1 = orient

        r = np.array([r.dot(x1),r.dot(y1),r.dot(z1)])

        result =0

        for i in range(0,self.M):
            for j in range(0,self.N):
                result += self.matrix[i][j].intensity(r)
        return result

    def LED_set(self, i,j, I=False, m=False, location=False, theta=False, phi=False):
        #False: remain default setting
        if I is not False:
            self.matrix[i][j].I = I
        if m is not False:
            self.matrix[i][j].m = m
        if location is not False:
            self.matrixp[i][j].location = location
        if theta is not False or phi is not False:
            t = 0 if theta is False else theta
            p = 0 if phi is False else phi
            self.matrix[i][j].set_orientation_angle(t,p)

    def LEDs_angle_reset(self):
        for i in range(0, self.M):
            for j in range(0, self.N):
                self.matrix[i][j].set_orientation_angle(0,0)

    def plot_intensity_data(self, xmin, xmax, ymin, ymax, z, xgrid=1E-3,ygrid=1E-3,position=np.array([0,0,0]), orient = np.array([0,0,1])):
        xi = np.arange(xmin,xmax+xgrid,xgrid) 
        yi = np.arange(ymin,ymax+ygrid,ygrid) 
        xdata,ydata=np.meshgrid(xi, yi)
        zdata = np.zeros(shape =xdata.shape)

        for i in xdata.shape[0]:
            for j in xdata.shape[1]:
                tar = np.array([xdata[i][j],ydata[i][j],z])
                zdata[i][j] = self.intensity(tar, position=position, orient =orient)
        return xdata, ydata ,zdata

#===================================================
# Morena et al, Sparrow criterion methods functions
# This methods does not use LED class only calculate x,y array for each LED for linear and rectangular array. 
#===================================================
def morena_coefficient(m,N,M=1,shape="L", approx = False):
# Varaible check 


#---------------------------------------------

    result = 0
# Calculation functions
    def morena_linear(x,N):
        result = 0
        for i in range(1,N+1):
            result += (1-(m+3)*(N+1-2*i)**2 * (x**2)/4)*((N+1-2*i)**2 * (x**2)/4 +1)**(-(m+6)/2)
        return result

    def morena_rectangular(x, N, M):
        result = 0
        for i in range(1, N+1):
            for k in range(1, M+1):
                result +=(1-((m+3)*(N+1-2*i)**2-(M+1-2*j)**2)*(x**2/4))*(((N+1-2*i)**2 + (M+1-2*j)**2)*(x**2/4)+1)**(-(m+6)/2)

# "L" : linear case
# "R" : Rectangular case

    if shape == "L":
        if approx == True and N>4 and m >30: #approximation coefficient, fast
            result = math.sqrt(3.2773/m+4.2539)
        else:
            if N%2 ==0: # even LED
                if N ==2: # n=2 case is simple
                    result = math.sqrt(4/(m+3))
                else:
                    r = op.root_scalar(lambda x : morena_linear(x,N), bracket=[0,1], method = "brentq")
                    result = r.root
            else :
                r = op.minimize_scalar(lambda x: morena_linear(x,N), bounds=(0,1), method="bounded")
                result = r.x
    elif shape == "R":
        if approx == True and (N>4 and M>4) and m >30:
            result = math.sqrt(1.2125/(m-3.349))
        else:
            if N == M and N==2:
                result = math.sqrt(4/(m+2))
            else:
                try :
                    r = op.root_scalar(lambda x: morena_rectangular(x,N,M), bracket=[0,1],method="brentq")
                    if r.converged == False:
                        r = op.minimize_scalar(lambda x: morena_rectangular(x,l,n,m), bounds=(0,1), method="bounded")
                        result = r.x
                    else:
                        result = r.root
                except:
                    r = op.minimize_scalar(lambda x: morena_rectangular(x,N,M), bounds=(0,1), method="bounded")
                    result = r.x
    return result


def morena_array(m,h,N,M=1,shape="L",approx=False,half=True):

# Varaible check 


#---------------------------------------------    

    if shape == "L":
        d = h * morena_coefficient(m,N,M,shape=shape, approx = approx)
        if half:
            if N%2 :
                array = np.fromfunction(lambda i,j: (j) *d, (1,int((N+1)/2)), dtype =int)
            else :
                array = np.fromfunction(lambda i,j : (1/2 +j)*d, (1,int(N/2)), dtype=int)
        else :
            array = np.fromfunction(lambda i,j: (-(N-1)/2+j)*d, (1,N), dtype=int)
        
        return [(N,1),d, (d*(N-1),0), array ]
    elif shape == "R":
        pass


#-------------------------------------------------------------


#===================================================
# Boundary enhancement method by Hyeon.
# This method does not use LED class only calculates x,y array for each LED for linear and rectangular array. 
#===================================================

def __intensity(x,t,phi, i, h ,m):
    if math.isclose(phi,0,rel_tol=np.finfo(float).eps):
        return h**m * i/(h**2 + (x-t)**2)**(m/2+1)
    else:
        sphi = math.copysign(1,phi)
        l = sphi*(x-h/math.tan(phi))
        sphi = math.copysign(1,phi)
        if sphi*t< l:
            messge = "Check phi, x, t\n sgn(phi) t > sgn(phi)*(x-h cot(phi))\n x:{}, t:{}, phi:{}".format(x,t,phi)
            raise ValueError(messge)
        r2 = (h**2 +(x-t)**2)**(m/2 +1)
        ang = (h * math.cos(phi) + (t-x)*math.sin(phi))**m

    return i *ang/r2 

def Ic(t,h,m,I, phi): 
        return (2* __intensity(t,0,phi,I,h,m))
def Ib(t,h,m,w,I,phi):
        return (__intensity(-t,w/2,-phi,I,h,m)+__intensity(t,w/2,phi,I,h,m))

def xe(h,w,m):
    r = op.root_scalar(lambda x:Ib(x,h,m,w,1, 0)/Ic(x,h,m,1, 0) -1, bracket=[0,w/2], method="brentq")
    return r.root

def xm(h,w,m,xe):
    L = Ib(w/2,h,m,w,1,0)-Ic(w/2,h,m,1,0)
    r= op.root_scalar(lambda x: Ib(x,h,m,w,1,0)-Ic(x,h,m,1,0)+L,bracket=[0, xe], method="brentq" )
    return r.root



def find_corres_BC(arr,xe,xm,h,w,m, append=True):
    def bc_find(x,xe,xm,h,w,m):
        if x<0 or x>w/2:
            raise ValueError("Argument 'x': must be in range [{},{}] \n x value= {}".format(0,w/2,x, x>0.0))
        elif math.isclose(x,xe,rel_tol=np.finfo(float).eps):
            return xe
        elif math.isclose(x,w/2, rel_tol=np.finfo(float).eps):
            return xm
        elif math.isclose(x,0,rel_tol=np.finfo(float).eps):
            return w/2
        elif xm < x and x < xe:
            # find the xc in [xs,w/2]'
            L = Ic(w/2,h,m,1,0)-Ib(w/2,h,m,w,1,0)
            r=op.root_scalar(lambda x: Ib(w/2,h,m,w,1,0)-Ic(w/2,h,m,1,0) -L, bracket=[xe,w/2], method="brentq")
            return r.root
        elif  x<xm:
            return w/2
        elif xe< x and x <w/2:
            # find the xc in [xm,xs]

            L = Ib(w/2,h,m,w,1,0)-Ic(w/2,h,m,1,0)
            r=op.root_scalar(lambda x: Ic(w/2,h,m,1,0)-Ib(w/2,h,m,w,1,0) -L, bracket=[xm,xe], method="brentq")
            return r.root
    
    correspoints = []
    for x in arr[0]:
        correspoints.append(bc_find(x,xe,xm,h,w,m))
    if append:
        return np.array([np.append(arr[0], np.array(correspoints))])
    else:
        return np.array([correspoints])


def H_matrix(m,w1,w2,h,I0=1):
    def morena_linear(m,h,w):
        n =2
        L  = morena_array(m,h,n,M=1,shape="L",approx=False,half=True)
        Lw = L
        while Lw[2]/2 < w/2: #need more condition. xm < x <xe
            L = Lw
            n += 1
            Lw = morena_array(m,h,n,M=1,shape="L",approx=False,half=True)
        return L


    xex = xe(h,w1,m)
    xey = xe(h,w1,m)
    xmx = xm(h,w2,m,xex)
    xmy = xm(h,w2,m,xey)

    xarr = morena_linear(m,h,w1)
    yarr = morena_linear(m,h,w2)

    xarray = find_corres_BC(xarr, xex,xmx,w1,m)
    yarray = find_corres_BC(yarr, xey,xmy,w2,m)

    return LEDmatrix(m,I0,xarray,yarray)
   


def generate_LEDmatrix(m,h,N,M=1, method="hyeon"):

    if method == "hyeon":
        pass
    elif method == "morenaL":
        pass
    elif method == "morenaR"

    return LEDmatrix