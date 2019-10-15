#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 12:46:33 2019

@author: lillian
"""
import numpy as np 
import matplotlib.pyplot as plt
import sympy as sy
from scipy.optimize import fsolve
import pandas as pd 
import matplotlib.pyplot as plt
import quantecon as qe
from scipy.optimize import brentq
from scipy.optimize import minimize

'''
Value Function Approach 
We try to solve for the maximizing function that approaches the steady state in optimal way (capital maximising outcome)
1. Set Grid: numers of future k that are taken into account to compute the value function 
2. Guess a solution 
3. Compute the indirect utility for each state 
4. Select maximizing vector of indirect utility as new Value Function 
5. Compare if Value Function improved
'''

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

def c(k,kk):
    return k**(1-theta)*1**theta+(1-delta)*k-kk

def u(c, h):
    return  np.log(c)-kappa*h
    #return  np.log(c) -kappa*((h**(1+(1/vue)))/(1+(1/vue)))
# =============================================================================
# 1.1 Value Function - bruce force iteration 
# =============================================================================

# 1. Step: Setting the grid 
N_k = 50
N_h=100
#N_h = 100 # h_grid = np.linspace(0.01, 1, N_h)
k_grid = np.linspace(0.01, 10, N_k)
h_grid = np.linspace(1, 1, N_h)

# 2. Step: Guess the Value function that we improve by iteration 
V = np.zeros((N_k))

### Empty policy functions
k_policy = np.zeros((N_k))
#c_policy = np.zeros((N_k))
c1_policy = np.zeros((N_k))
V_next = np.zeros((N_k))
m = np.zeros((N_k,N_k))
a = np.zeros((N_k,N_k))
negativ=-10000
i = 1 
j=0

'''Loop for: 
    1. Check if consumption>0  (non-negativity constraint)
    2. Construct indirect utility at k_i for each k_j
    3. Based on V construct V_next 
    4. Run Policy functions
    5. Check if V_next is a better guess than V, by sufficient iterations they become equal
    '''
while i<500:
    for ik, k in enumerate(k_grid): #ik is index number for the rows running from 0-->k_grid, while k in k_grid is an actual value
        for jk, kk in enumerate(k_grid): #jk is the index number for the columns
            c1=(y_func(k,h)+(1-delta)*k)-k_grid
            if c1[jk]>0:                        #checking consumption >0
                m[ik,:] = u(c[jk],1)+beta*V[jk]
            else:
                m[ik,jk] = negativ  
        V_next[ik] = np.nanmax(m[ik,:]) #value
        k_policy[ik] = np.argmax(m[ik,:])   #position
        #c_policy[ik] =  y_func(k,1) +(1-delta)*k - k_grid(k_policy[ik])  
        c1_policy[ik] =  y_func(k,1) +(1-delta)*k - k_policy[ik] 
    if abs(np.nanmax(V_next))-abs(np.nanmax(V))<0.001:  #  if np.allclose(V_next,V)==True:
        for i in range (0,49):
            if abs(V[i]-V_next[i])<0.05:
                j+=j
                if j>=45: 
                    print ('V-max found'+V_next)
                    break
    else:
        V_next=V
        i+=i
        
#Plotting graphs 
fig, ax = plt.subplots()
ax.plot(k_grid,V_next, label='V(k)')
plt.show()

### Plots
fig, ax = plt.subplots()
ax.plot(k_grid,V, label='V(k)')
plt.show()

fig, ax = plt.subplots()
ax.plot(k_grid,c1_policy, label='c(k)')
ax.plot(k_grid,k_policy, label='k1(k)')
ax.legend()
plt.show()
        
    
# =============================================================================
# 1.1 Value Function - bruce force iteration - !running!
# =============================================================================      
V_guess=np.zeros((N_k))
def bellman_operator(V,return_policies=False):
    V_next = np.zeros((N_k))
    m = np.zeros((N_k,N_k)) #why define again?
    
    for ik, k in enumerate(k_grid):
        for igh, gh in enumerate(h_grid):
            ### Fill up return matrix: All possible choices of k' (k_grid) and labor (gh) per each state (k)
            m[ik,:] =  u((y_func(k,igh) +(1-delta)*k - k_grid),gh) +beta*V
           
        V_next[ik] = np.nanmax(m[ik,:])     #wert   
        k_policy[ik] = np.argmax(m[ik,:])   #position
        c_policy[ik] =  y_func(k,1) +(1-delta)*k - k_policy[ik]         
    if return_policies==True:
        return V_next, k_policy, c_policy, m
    else:
        return V_next
#Compute fixed point in bellmann equation
qe.tic()
V = qe.compute_fixed_point(bellman_operator, V_guess, max_iter=1000, error_tol=0.001, print_skip=20)
V, g_k, g_c, m = bellman_operator(V, return_policies=True)
qe.toc()


### Plots
fig, ax = plt.subplots()
ax.plot(k_grid,V, label='V(k)')
plt.show()

fig, ax = plt.subplots()
ax.plot(k_grid,g_c, label='c(k)')
ax.plot(k_grid,g_k, label='k1(k)')
ax.legend()
plt.show()


# =============================================================================
# 1.1 Value Function - Monotonicity 
# =============================================================================
while i<500:
    for ik, k in enumerate(k_grid): #ik is index number for the rows running from 0-->k_grid, while k in k_grid is an actual value
        for jk, kk in enumerate(k_grid): #jk is the index number for the columns
            c1=(y_func(k,h)+(1-delta)*k)-k_grid
            if c1[jk]>0:                        #checking consumption >0
                m[ik,:] = u(c[jk],1)+beta*V[jk]
            else:
                m[ik,jk] = negativ  
            m[ik,:]=np.where(jk<ik, None, m[ik,:]) #deleting all values that are non-relevant for the policy function
        V_next[ik] = np.nanmax(m[ik,:]) #value
        k_policy[ik] = np.argmax(m[ik,:])   #position
        #c_policy[ik] =  y_func(k,1) +(1-delta)*k - k_grid(k_policy[ik])  
        c1_policy[ik] =  y_func(k,1) +(1-delta)*k - k_policy[ik] 
    if abs(np.nanmax(V_next))-abs(np.nanmax(V))<0.001:  #  if np.allclose(V_next,V)==True:
        for i in range (0,49):
            if abs(V[i]-V_next[i])<0.05:
                j+=j
                if j>=45: 
                    print ('V-max found'+V_next)
                    break
    else:
        V_next=V
        i+=i