import numpy as np


mp =1 #mass of particle 
   
def Force(x,k,V,dH):                   # Force at position x and surface k, considering the potential and derivative of the hamiltonian
     eigval,eigvect = np.linalg.eigh(V(x)) # finding eig funcions of V 
     f = - eigvect[:,k].transpose()@ dH(x) @ eigvect[:,k]
     return np.asscalar(f)

def accel(x,k,V,H):                    # acceleration of particle at position x and surface k 
    a = Force(x,k,V,H)/mp
    return a
  
def verlet_pos(xi,vi,dt,k,V,H):        # finds new posion xf from initial position xi, intial veloxity vi, on surface k during a time step of length dt
        xf = xi + vi *dt+accel(xi,k,V,H)*dt**2*0.5
        return xf
def verlet_vel(xi,xf,vi,dt,k,V,H):    # finds new velocity vf from intial position xi, the finial position xf, the initial velocity vi, on surface k, time step dt 
        vf = vi+0.5*dt*(accel(xi,k,V,H)+accel(xf,k,V,H))
        return vf



