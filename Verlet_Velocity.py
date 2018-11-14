import numpy as np
import math

   
def Force(x,surface,V,H):
     eigval,eigvect = np.linalg.eigh(V(x))
     f = - eigvect[:,surface].transpose()@ H(x) @ eigvect[:,surface]
     return np.asscalar(f)

def accel(x,surface,V,H):
    a = Force(x,surface,V,H)/2000
    return a
  
def verlet_pos(xi,vi,dt,k,V,H):
        xf = xi + vi *dt+accel(xi,k,V,H)*dt**2*0.5
        return xf
def verlet_vel(xi,xf,vi,dt,k,V,H):
        vf = vi+0.5*dt*(accel(xi,k,V,H)+accel(xf,k,V,H))
        return vf



