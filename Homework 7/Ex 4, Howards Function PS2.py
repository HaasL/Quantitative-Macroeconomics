#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 09:16:08 2019

@author: lillian
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 10:31:00 2019

@author: lillian
"""

import numpy as np 
import matplotlib.pyplot as plt
import quantecon as qe
import time

# =============================================================================
# Excercise 1 - Set up
# =============================================================================
# Functions

def u(c):
    if theta==1: 
        u=np.log(c)
    else: 
        u=c**(1-theta)/(1-theta)
    return u

def func_extrapol(x1,x2,y1,y2,x):
    m = (y2-y1)/(x2-x1)
    return y1 + m*(x-x1)

def makegrid(x1,x2,nc,curv):
    scale =x2-x1
    grd[0]=x1
    grd[nc-1]=x2
    #a=nc-1
    for i in range (1,nc):
        grd[i]=x1+scale*((i-1.0)/(nc-1.0))**curv
    return grd

def margutil(c):
    c=max(c,min_cons)
    muc = c**(-theta)
    return muc


def func_intp(x,func,xp):
    #chkout=false
    n = len(x)
    if xp>x(n):
        fv=func_extrapol(x[n-1],x[n],func[n-1],func[n],xp)
    elif xp<x[1]:
        fv=func_extrapol(x[1],x[2],func[1],func[2],xp)
    #else:
        #fv = interp1(x,func,xp)
    return fv


def evalvp(cons,x): 
    vpp1 = np.zeros(ne,1)
    for ec in range (1,ne):
        xp1 = (x-cons)*(1+r)/(1+g)+epsi[ec]
        vpfun = margutil(c) #check if this is the vpfun necessary 
        vpp1[ec] = func_intp[a_grid,vpfun,xp1]
        vpp1[ec] = vpp1[ec]*probepsi[ec] #dot product
    return sum(vpp1)
        

#Settings 
eps=0.000001
theta = 1
rho=0.03
beta = 1/(1+rho)
r = 0.02
g=0.01 
kappa = 5.24
vue = 2

maxit = 100 
tol=1e-4
nt=1100    
dt=100     
min_cons=1.0e-08


#Grid Settings 
N=99
N_a = 100           #x grid points 
amax = 30      
amin = np.sqrt(eps)
curv=3.0            #curvature of grid
grd=np.zeros(N_a)
a_grid=makegrid(amin, amax, N_a, curv) #check for code  or/and transpose


#Income shocks
ne = 7
varepsi = 0.01
muepsi = -varepsi/2
[epsi,probepsi] = qe.quad.qnwnorm(ne,muepsi,varepsi)
epsi=np.exp(epsi)
#if abs(sum(np.dot(epsi,probepsi))-1.0)>np.sqrt(eps):
#    print('random numbers fucked up')
    


# =============================================================================
# Excercise 4 - Brute force method 
# =============================================================================
negativ=-10000
i = 1 
j=0

# 2. Step: Guess the Value function     
V = np.zeros((N_a,ne))
V_guess=np.zeros ((N_a,ne)) #guess
V_new=np.zeros ((N_a,ne))
V_result=np.empty ((N_a,ne))

# Empty policy functions
a_policy = np.zeros((N_a,ne))
a_newpolicy = np.zeros((N_a,ne))
c_policy = np.zeros((N_a,ne))
V_next = np.zeros((N_a,ne))
m = np.zeros((N_a,N_a,ne))
m1 = np.zeros((N_a,N_a))
a = np.zeros((N_a,N_a,ne))
Chi = np.zeros((N_a,N_a,ne)) 
Chi_s = np.zeros((N_a,N_a))    


#3. Step: Creating the m-matrix
start_time = time.time() 
for ia, a in enumerate(a_grid): 
    for ja, aa in enumerate(a_grid): 
        for ka, aaa in enumerate(epsi):
            c=(aa/(1+r))+aaa +a 
            if aa>=0: #Check for budget constraint 
                m[ia,ja,ka] = u(c)      #checking consumption >0
            else:
                m[ia,ja,ka] = negativ         
              
while i<1: 
    if i<=1: 
        V=V_guess 
    else: 
        V=V_next
        V_next=np.zeros((N_a,ne))
    for ia, a in enumerate(a_grid):
        for ka, aaa in enumerate(epsi):
            for ja, aa in enumerate(a_grid):
                Chi[ia,ja,ka]=m[ia,ja,ka]+beta*(1+g)**(1-theta)(1+r)*V[ia,ka]
            V_new[ia,ka] = np.nanmax(Chi[ia,ka,:])    
            V_next[ia,ka]= V[ia,ka]*probepsi[ka]
            #b=np.mean(sum(V_next[ia,:]))
            V_next[ia,ka] =  np.mean(sum(V_next[ia,:])) #value
            a_policy[ia,ka] = np.argmax(Chi[ia,ka,:]) #position Policy Function
    if np.allclose(V_next,V): 
        V_result=V_new
        break 
    else:        
        i+=1   


for ia, a in enumerate(a_grid):    
     for ka, aaa in enumerate(epsi):
         iaa=ia+1
         if iaa<N_a:
             c_policy[ia,ka]=a_policy[iaa,ka]/(1+r)+epsi[ka]+a_policy[ia,ka]
         #if iaa==N_a:
             #c_policy[ia,ka]=epsi[ka]+a_policy[ia,ka]
V_use=V_next             
a=5
#4. Generating a Value Function conditional on last policy function 
Runs=[a,a,a,a,a,a] #Number of iterations for new Value functions before a new Policy Function is generated
for A in Runs:
    j+=1
    for i in range (0,A):
        
        for ia, a in enumerate(a_grid): #Creating a new Value Function 
            for ka, aaa in enumerate(epsi):
                aa= a_policy[ia,ka]
                aa=int(aa)
                iaa=ia+1
                if iaa< 100: 
                    if a_grid[iaa]>=0: 
                        aaa=a_grid[iaa]
                        u=np.log(aa/(1+r)+aaa +a) 
                        V_next[ia,ka] = u +beta*(1+g)**(1-theta)*(1+r)*V_use[iaa,ka]  #checking consumption >0
                    else:
                        V_next[ia,ka] = negativ
                else:  
                    aaa=0
                    u=np.log(aa/(1+r)+aaa +a) 
                    V_next[ia,ka] = u 
                    
                  
                
        if np.allclose(V_next,V_use):  #Checking if convergence was achieved
            V_result=V_next
            print('Optimal Policy was computed in', A, i)
            break
        else: 
            V_use=V_next
    #After 'A' iterations: Create new policy function (capital) 
    for ia, a in enumerate(a_grid):
        for ka, aaa in enumerate(epsi):
            for ja, aa in enumerate(a_grid):
                Chi[ia,ja,ka]=m[ia,ja,ka]+beta*(1+g)**(1-theta)*(1+r)*V[ia,ka]
            V_new[ia,ka] = np.nanmax(Chi[ia,ka,:])    
            V_next[ia,ka]= V[ia,ka]*probepsi[ka]
            #b=np.mean(sum(V_next[ia,:]))
            V_next[ia,ka] =  np.mean(sum(V_next[ia,:])) #value
            a_newpolicy[ia,ka] = np.argmax(Chi[ia,ka,:]) #position Policy Function
    if np.allclose(V_next,V): 
        V_result=V_new
        print('Optimal Policy was computed in Policy Check', A, ia)
        break 
    elif (a_newpolicy==a_policy).all():   #Checking if convergence was achieved
        print('Optimal Policy was computed along new Policy')
        break
    else: 
        a_newpolicy=a_policy           
        
print(j,'iterations are necessary to converge')
V_result=V_next    



             
        
        
end_time = time.time()     
print("total time taken this loop: ", end_time - start_time)   
print('total number of iterations: ', i)     
        
#Plotting graphs 
fig, ax = plt.subplots()
ax.plot(a_grid,V_result, label='V(k)')
plt.title ('Value Function Iteration')
plt.legend()
plt.xlabel('Time') 
plt.show()

### Plots
fig, ax = plt.subplots()
ax.plot(a_grid,c_policy, label='c(k)')
ax.plot(a_grid,a_policy, label='a(k)')
plt.title ('Policy Functions')
plt.legend()
plt.xlabel('Time') 
plt.show()








                

