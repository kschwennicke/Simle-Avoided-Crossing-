import numpy as np


def C_dot_k(C,x,v,k,V,d):
    eigvals = np.linalg.eigvalsh(V(x))
    V = np.matrix([[eigvals[0],0],[0,eigvals[1]]])
    a0 = C[0]*(V[k,0]-complex(0,v *d(x)[k,0]))
    a1=C[1]*(V[k,1]-complex(0,v *d(x)[k,1]))
    a = (a0+a1)*complex(0,-1)
    return a

def Ck (C,xi,x1_2,xf,vi,v1_2,vf,k,dt,V,d):
        C0 = C[0]
        C1 =C[1]
        if k ==0:
            
            k1 = dt*C_dot_k(C,xi,vi,k,V,d)
            C1_2 = [C0+k1/2,C1]
            k2 = dt*C_dot_k(C1_2,x1_2,v1_2,k,V,d)
            C1_2 = [C0+k2/2,C1]
            k3 = dt*C_dot_k(C1_2,x1_2,v1_2,k,V,d)
            Cf = [C0+k3,C1]
            k4 = dt*C_dot_k(C1_2,xf,vf,k,V,d)

            Ck = C0 + 1/6*(k1+2*k2+2*k3+k4)
        
        
        if k ==1:
            
            k1 = dt*C_dot_k(C,xi,vi,k,V,d)
            C1_2 = [C0,C1+k1/2]
            k2 = dt*C_dot_k(C1_2,x1_2,v1_2,k,V,d)
            C1_2 = [C0,C1+k2/2]
            k3 = dt*C_dot_k(C1_2,x1_2,v1_2,k,V,d)
            Cf = [C0,k3+C1]
            k4 = dt*C_dot_k(C1_2,xf,vf,k,V,d)

            Ck = C1 + 1/6*(k1+2*k2+2*k3+k4)
   
        return Ck



    


        
        
        

