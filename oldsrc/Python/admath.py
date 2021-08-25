import numpy as np
import scipy as scp
from scipy import special as sp
from scipy import optimize as op
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
import random
import itertools

d= 20E-3
m = 80
width = 200E-3
n =400

intensity = 100

def kernel(x,t,d,m):
    return d**m /(d**2 + (x-t)**2)**(m/2 +1)

def fg(i,j):
    return width * pair[1][j] * kernel(width/2 *pair[0][i], width/2 *pair[0][j],d,m)

def f(i,j):
    return width * kernel(width/2 *pair[0][i], width/2 *pair[0][j],d,m)

dx = width/(n+1)
pair = sp.roots_legendre(n)
Fg = np.fromfunction(fg, (n,n), dtype=int)
F = np.fromfunction(f, (n,n), dtype=int)
Intense = intensity * np.ones(n)
solg = op.nnls(F,Intense)[0] 
sol2 =op.lsq_linear(F,Intense, bounds=(0,1),lsmr_tol="auto",verbose=1).x

plt.plot(width/2 *pair[0], sol)
plt.axis([-width/2, width/2, 0, 40])
plt.show()
plt.plot(width/2 *pair[0], sol2)
plt.axis([-width/2, width/2, 0, 1.2])
plt.show()

#-------------------------------------------------------------------------------
# Genetic algorithm

@dataclass(order=True)
class Obj:
    bvec: field(default_factory= np.ndarray)
    eval: np.double
    length: int
    probability: np.double



class binary_set:
    def __init__(self,n,m):
        self.soldim=n
        self.spacedim =m
        self.obj =self.generate()
        
    def generate(self):
        obj = list()
        for i in range(0,self.spacedim):
            a=np.empty([1,self.soldim], dtype=np.uint64)
            b=0
            while b == 0:
                c = np.random.randint(2,size=self.soldim)
                for j in range(0,i):
                    b += np.array_equal(c, obj[i].bvec)
                if b == 0  :
                    a= c.copy()
                else :
                    b = 0
            obj.append(Obj(a,0.0,self.soldim, 0.0))
        return  obj

class Binary_solution:
    def __init__(self, mp, mu, T, width, d, m, n, k):
        #check m > g(k) <- possible next generation
        self.obj = binary_set(n,m)
        self.F = np.fromfunction(f, (2*n,2*n), dtype=np.uint64)
        self.soldim=n
        self.spacedim = m
        self.picknum =k

        #Probability-----------------
        self.T = T #Boltzmann temperature
        self.mu = mu #reducing coefficient
        self.mp = mp # mutation probaility
        #---------------------

    def ordering(self):
        self.obj = sorted(self.obj, key=lambda x: x.eval)
        
    def set_eval(self,obj):
        v = np.concatenate((self.obj.bvec, np.flip(self.obj.bvec)))
        self.obj.eval = np.std(np.dot(F,v))

    def eval_all(self):
        for i in range(0,self.spacedim):
            self.set_eval(self.obj[i])

    def set_probaility(self):
        b=0
        for i in range(0,self.spacedim):
            b += np.exp(self.obj[i].eval/(sp.k*self.T))
        for i in range(0,self.spacedim):
            self.obj[i].probability = np.exp(self.obj[i].eval/(sp.k*self.T))/b

    def reduce_temperature(self):
        self.T = self.T / self.mu
    
    def select_parent(self): #Boltzmann Selection
        self.parent = list()
        for i in range(0,self.picknum):
            p = random.random()
            for j in range(0, self.spacedim):
                p = p - self.obj[j].probability
                if p <0 :
                    self.parent.append(self.obj[j])
    
    def formming_next_generation(self):

        self.obj = list()
        comset = list(itertools.combinations(self.parent,2))

        for i in range(0,self.picknum):
            self.obj.append(self.parent[i])
            self.obj.append(self.mutation(self.parent[i]))

        for com in comset:

            #one Crossover
            cros = self.crossover(0,com[0], com[1])
            self.obj.append(cros[0])
            self.obj.append(cros[1])
            self.obj.append(self.mutation(cros[0]))
            self.obj.append(self.mutation(cros[1]))
            #two Crossover
            cros = self.crossover(1,com[0], com[1])
            self.obj.append(cros[0])
            self.obj.append(cros[1])
            self.obj.append(self.mutation(cros[0]))
            self.obj.append(self.mutation(cros[1]))
            #Uniform Crossover
            cros = self.crossover(2,com[0], com[1])
            self.obj.append(cros[0])
            self.obj.append(cros[1])
            self.obj.append(self.mutation(cros[0]))
            self.obj.append(self.mutation(cros[1]))

    def crossover(self, case, obj1, obj2): #return 2 obj, 0: one, 1: two, 2: uniform
        if case ==0: #0: one
            l = self.soldim % 2
            if l ==0: # n = even

            else:                  # n = odd

        elif case == 1: #1: two
            l = self.soldim % 3
            if l == 0:

            elif l ==1:

            elif l ==2:

        elif case == 2: #2: uniform
            l = random.randint(3, self.soldim)
            flag = np.sort(np.random.choice(self.soldim-1, l, replace=False)+1)


    def mutation(self, obj):
        p = 0
        if self.mp ==0:
            p = 1/self.soldim
        else:
            p =self.p
        a = obj.bvec.copy()
        b= np.zeros(self.soldim, dtype=np.uint64)
        for i in range(0,self.soldim):
            if random.random() < p :
                b[i] = 1
                
        return Obj(np.logical_xor(a,b).astype(np.uint64),0.0,self.soldim, 0.0)

