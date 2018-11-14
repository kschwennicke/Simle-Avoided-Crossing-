import math
import numpy as np
import matplotlib.pyplot as plt

####Plotting Adiabatic Potential Funtion####

A = 0.01
B = 1.6
C= 0.005
D =1.0

def V11(x):
    if x >0:
        a = A*(1.- math.exp(-B*x))

    if x <= 0:
        a = - A*(1.-math.exp(B*x))
    if x == 0:
        a =0
    return a

def V22(x):
    b = -V11(x)
    return b

def V12(x):
    c= C*math.exp(-D*x**2)
    return c


def Vmatrix (x):
    V = np.matrix([[V11(x), V12(x)], [V12(x), V22(x)]])
    return V
def diag(m):
    w = np.linalg.eigvalsh(m)
    matrix = np.matrix([[w[0],0],[0,w[1]]])
    return matrix

def dH11 (x):
    if x>0:
        y = A*B*np.exp(-B*x)

    if x<0:
    
        y = A*B*np.exp (B*x)
    if x ==0:
        y = A*B

    return y

def dH22(x):
    y = - dH11(x)
    return y
def dH12 (x):
    y = - C*2*D*x *np.exp(-D*x**2)
    return y

def dHmatrix(x):
    y = np.matrix([[dH11(x),dH12(x)],[dH12(x),dH22(x)]])
    return y

def dmatrix (x):
    val,vect = np.linalg.eigh(Vmatrix(x))
    y  = vect[:,0].transpose()@ dHmatrix(x)@ vect[:,1]
    w =np.asscalar(y)/(val[1]-val[0])
    z = np.matrix([[0,w],[-w,0]])
    return z

def Force(x,surface):
     eigval,eigvect = np.linalg.eigh(Vmatrix(x))
     f = - eigvect[:,surface].transpose()@ dHmatrix(x) @ eigvect[:,surface]
     return np.asscalar(f)

def accel(x,surface,mp):
    a = Force(x,surface)/mp
    return a
  


    


