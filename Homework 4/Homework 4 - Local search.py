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
# 1.4 Value Function - Local Search
# =============================================================================
# 1. Step: Setting the grid 
N_k = 100

k_grid = np.linspace(0.5, 70, N_k)


# 2. Step: Guess the Value function that we improve by iteration 
V = np.zeros((N_k))
V_guess=np.zeros ((N_k))
V_result=np.empty ((N_k))
V_next = np.zeros((N_k))

### Empty policy functions
k_policy = np.zeros((N_k))
c_policy = np.zeros((N_k))
V_next = np.zeros((N_k))
m = np.zeros((N_k,N_k))
Chi = np.zeros((N_k,N_k))
negativ=-10000
k_local=0
i = 1 
error=102 #102 is the 'right' search radius around k_t

#3. Step: Creating the m-Matrix
for ik, k in enumerate(k_grid): #ik is index number for the rows running from 0-->k_grid, while k in k_grid is an actual value
    for jk, kk in enumerate(k_grid): #jk is the index number for the columns
        if kk<y_func(k,h)+(1-delta)*k:
            m[ik,jk] = u(y_func(k,h) +(1-delta)*k - kk) #checking consumption >0
        else:
            m[ik,jk] = negativ


#4. Step: Finding the Value Function
while i<10000: 
    start_time = time.time()
    if i<=1: 
        V=V_guess 
    else: 
        V=V_next
        V_next=np.zeros((N_k))
    for ik, k in enumerate(k_grid): #ik is index number for the rows running from 0-->k_grid, while k in k_grid is an actual value
        if ik<1:
            for jk, kk in enumerate(k_grid): #jk is the index number for the columns
                Chi[ik,jk]=m[ik,jk]+beta*V[jk] #checking consumption >0
        else:                   #Creating for all ik>=1 Chi-Values in a certain radius(=error)
            low=k_local-error
            high=k_local+error
            if low<0:   #Controlling for k-values in k_grid
                low=0 
            if high>70: 
                high=70
            for jk in range (low, high): #radius or only considering right side due to monotonicity
                Chi[ik,jk]=m[ik,jk]+beta*V[jk]        
        V_next[ik] = np.nanmax(Chi[ik,:]) #value
        k_policy[ik] = np.argmax(Chi[ik,:])   #position
        c_policy[ik] = y_func(k,1) +(1-delta)*k - k_policy[ik]   
        k_local=k_policy[ik]   
    if np.allclose(V_next,V): 
        V_result=V_next
        break 
    else:
        i+=1
end_time = time.time()  

print("total time taken this loop: ", end_time - start_time)   
print('total number of iterations: ', i)   

#Plotting Value Function
fig, ax = plt.subplots()
ax.plot(k_grid,V_result, label='V(k)')
plt.title ('Value Function Iteration')
plt.legend()
plt.xlabel('Time') 
plt.show()

### Ploting Policy Functions
fig, ax = plt.subplots()
ax.plot(k_grid,c_policy, label='c(k)')
ax.plot(k_grid,k_policy, label='k(k)')
plt.title ('Policy Functions')
plt.legend()
plt.xlabel('Time') 
plt.show()