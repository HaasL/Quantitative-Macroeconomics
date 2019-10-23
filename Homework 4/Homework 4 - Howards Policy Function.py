#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:34:10 2019

@author: lillian
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 18:36:01 2019

@author: lillian
"""

import numpy as np 
import matplotlib.pyplot as plt
import time

# =============================================================================
# Excercise 1 - Set up
# =============================================================================
#Given data by PS4 
theta = 0.679
beta = 0.988
delta = 0.013
h=1 
kappa = 5.24
vue = 2

# Functions
def y_func(k, h):
    return k**(1-theta)*h**theta

def c(k,kk,h):
    return k**(1-theta)*h**theta+(1-delta)*k-kk

def u(c):
    return  np.log(c)
# =============================================================================
# 1.4 Value Function - Howards Policy Function
# =============================================================================

# 1. Step: Setting the grid 
N_k = 200

k_grid = np.linspace(0.5, 70, N_k)

# 2. Step: Guess the Value function that we improve by iteration 
V = np.zeros((N_k))
V_guess=np.zeros ((N_k))
V_result=np.empty ((N_k))
V_next = np.zeros((N_k))
V_use= np.zeros((N_k))

### Empty policy functions
k_policy = np.zeros((N_k))
k_newpolicy = np.zeros((N_k))
c_policy = np.zeros((N_k))
m = np.zeros((N_k,N_k))
Chi = np.zeros((N_k,N_k))
negativ=-10000

i = 1 
j=0


#3. Creating the m matrix for all possible capital options
for ik, k in enumerate(k_grid): 
    for jk, kk in enumerate(k_grid): 
        if kk<y_func(k,h)+(1-delta)*k:
            m[ik,jk] = u(y_func(k,h) +(1-delta)*k - kk)    #checking consumption >0
        else:
            m[ik,jk] = negativ

#Iterate for a 'i' times the Value Function
#This is done to create starting values for the policy function iteration            
while i<100: 
    start_time = time.time()
    if i<=1: 
        V=V_guess 
    else: 
        V=V_next
        V_next=np.zeros((N_k))
    for ik, k in enumerate(k_grid): 
        for jk, kk in enumerate(k_grid): 
                Chi[ik,jk]=m[ik,jk]+beta*V[jk]
        V_next[ik] = np.nanmax(Chi[ik,:]) #value
        k_policy[ik] = np.argmax(Chi[ik,:])   #position Policy Function
    if np.allclose(V_next,V): 
        V_result=V_next
        break 
    else:
        i+=1

V_use=V_next        
#4. Generating a Value Function conditional on last policy function 

Runs=[50,20,10,5] #Number of iterations for new Value functions before a new Policy Function is generated
for A in Runs:
    for i in range (0,A):
        for ik, k in enumerate(k_grid): #Creating a new Value Function 
            kk=k_policy[ik]
            kk=int(kk)
            if k_grid[kk]<y_func(k,h)+(1-delta)*k:
                kkk=k_grid[kk]
                u=np.log(k**(1-theta)*h**theta+(1-delta)*k-kkk) 
                V_next[ik] = u +beta*V_use[kk]   #checking consumption >0
            else:
                V_next[ik] = negativ
        if np.allclose(V_next,V_use):  #Checking if convergence was achieved
            V_result=V_next
            print('Optimal Policy was computed in', A, i)
            break
        else: 
            V_use=V_next
    #After 'A' iterations: Create new policy function (capital) 
    for ik, k in enumerate(k_grid): 
        for jk, kk in enumerate(k_grid): 
                Chi[ik,jk]=m[ik,jk]+beta*V_use[jk]
        k_newpolicy[ik] = np.argmax(Chi[ik,:]) #position
        V_next[ik] = np.nanmax(Chi[ik,:]) #value
        c_policy[ik] =  y_func(k,h) +(1-delta)*k - k_newpolicy[ik]    
    if np.allclose(V_next,V_use):        #Checking if convergence was achieved
        V_result=V_next
        print('Optimal Policy was computed in Policy Check', A, ik)
        break
    else:
        V_use=V_next
    if (k_newpolicy==k_policy).all():    #Checking if convergence was achieved
        print('Optimal Policy was computed along new Policy')
        break   
    else: 
        k_newpolicy=k_policy             #Checking if convergence was achieved
        if A==5: 
            print ('More iteration are necessary')
    print (A)

V_result=V_next    
    
end_time = time.time()    
print("total time taken this loop: ", end_time - start_time)    
print('total number of iterations: ', i)     
        
#Plotting Value Functions
fig, ax = plt.subplots()
ax.plot(k_grid,V_result, label='V(k)')
plt.title ('Value Function Iteration')
plt.legend()
plt.xlabel('Time') 
plt.show()

### Ploting Policy Function
fig, ax = plt.subplots()
ax.plot(k_grid,c_policy, label='c(k)')
ax.plot(k_grid,k_policy, label='k(k)')
plt.title ('Policy Functions')
plt.legend()
plt.xlabel('Time') 
plt.show()
