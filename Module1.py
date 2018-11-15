import numpy as np
import RK4                                                  # importing RK4 method and diff equaion
import InitPos as guass                                     # import guassian disrubition 
import Simple_avoided_crossing as PE                        # importing our potential surface 
import Verlet_Velocity as vv                                # importing out Verlet Velocity module 
import matplotlib.pyplot as plt


                                                    
dt = 1                                                      # time step
alpha = 0.25; sigr = np.sqrt(1/(2*alpha)); r0 =-7; n = 10   # gauss distrubution 
g = 0                                                       # ground state
e =1                                                        # excited state
V = PE.Vmatrix                                              # potential matrix 
dH = PE.dHmatrix                                            # derivative of potential matrix 
d = PE.dmatrix                                               #coupling matrix 

P0av =[]                                                  # "average population" of ground state (not entirely correct)
P1av =[]                                                  # average poulation" of excited state (not entirely correct)
for j in range (n):                                       # loop for each trajectory  
     x = []
     v = []
     a = []
     t = []
     P0 =[]
     P1 =[]
     norm =[]
# itial conditions
     Co = np.array([1,0])
     
     x.append(np.random.choice(guass.ipos(r0,sigr,n)))     
     v.append(.1)
     P0.append(1)
     P1.append(0)
     k = g                                                      
     a.append(vv.accel(x[0],k,V,dH))
     t.append(0)
    

#loop will be replaced/scrapped later
     
     for m in range(1,3000):
          xi = x[m-1] 
          vi = v[m-1]
          ti = t[m-1]

# final conditiona
          tf = t[m-1]+dt
          xf =vv.verlet_pos(xi,vi,dt,k,V,dH)
          vf =vv.verlet_vel(xi,xf,vi,dt,k,V,dH)
          x1_2 = vv.verlet_pos(xi,vi,dt/2,k,V,dH)
          v1_2 =vv.verlet_vel(xi,x1_2,vi,dt/2,k,V,dH)
          C0 = RK4.Ck(Co,xi,x1_2,xf,vi,v1_2,vf,0,dt,V,d)
          C1 = RK4.Ck(Co,xi,x1_2,xf,vi,v1_2,vf,1,dt,V,d)
          Co = np.array([C0,C1])
          w = np.linalg.eigvalsh(V(xf))
          U = np.matrix.([[w[0],0],[0,w[1]])
          p0 = U[0,0]*Co[0]*np.conj(U[0,0]*Co[0])+U[0,1]*Co[1]*np.conj(U[0,1]*Co[1])
          p1 = U[1,0]*Co[0]*np.conj(U[0,0]*Co[0])+U[0,1]*Co[1]*np.conj(U[1,1]*Co[1])
          norm.append( p0+p1)
          P0.append(p0)
          P1.append(p1)
          x.append(xf)
          v.append(vf)
          t.append(tf)

     # here I have v and d so I can calculate electronic dynamics and then hoping
     #     RK4.RK4()
     
     
          #if x[m] > 10:
              # break
          #elif x[m] < -10:
              # break
     #plt.figure(1)
     #plt.subplot(212)
     #plt.plot(x,t)
     
     
     
     if j == 0:
          P0av = P0
          P1av = P1
     else:
          for i in range(len(P0)):
               P0av[i] = P0av[i]+P0[i]
               P1av[i]=P1av[i]+P1[i]
for i in range(len(P0av)):
     P0av[i]=  P0av[i]/n
     P1av[i]= P1av[i]/n

#plt.figure(1)
#plt.subplot(211)
plt.plot(t,P0av,label='P0')
plt.plot(t,P1av,label ='P1')
#plt.plot(t,norm)
#plt.legend()
plt.ylabel('Probability')
    
     





plt.show()
    



