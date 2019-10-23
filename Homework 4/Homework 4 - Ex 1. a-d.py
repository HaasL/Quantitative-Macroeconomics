#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 19:38:37 2019

@author: lillian
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 12:46:33 2019

@author: lillian
"""
import numpy as np 
import matplotlib.pyplot as plt
import quantecon as qe
import time


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

def c(k,kk,h):
    return k**(1-theta)*h**theta+(1-delta)*k-kk

def u(c):
    return  np.log(c)
# =============================================================================
# 1.1 Value Function - bruce force iteration 
# =============================================================================
'''Loop for: 
    1. Check if consumption>0  (non-negativity constraint)
    2. Construct indirect utility at k_i for each k_j
    3. Based on V construct V_next 
    4. Run Policy functions
    5. Check if V_next is a better guess than V, by sufficient iterations they become equal
'''

# 1. Step: Setting the grid 
N_k = 100
k_grid = np.linspace(0.5, 70, N_k)

# 2. Step: Guess the Value function that we improve by iteration 
V = np.zeros((N_k))
V_guess=np.zeros ((N_k))
V_result=np.empty ((N_k))

### Empty policy functions
k_policy = np.zeros((N_k))
c_policy = np.zeros((N_k))
c_policy = np.zeros((N_k))
V_next = np.zeros((N_k))
m = np.zeros((N_k,N_k))
a = np.zeros((N_k,N_k))
Chi = np.zeros((N_k,N_k))

negativ=-10000
i = 1 
j=0

#3. Step: Creating the m-matrix
for ik, k in enumerate(k_grid): 
    for jk, kk in enumerate(k_grid): 
        if kk<y_func(k,h)+(1-delta)*k:
            m[ik,jk] = u(y_func(k,h) +(1-delta)*k - kk)    #checking consumption >0
        else:
            m[ik,jk] = negativ

#4. Step: Value Function iteration
while i<10000: 
    qe.tic()
    start_time = time.time()
    if i<=1: 
        V=V_guess 
    else: 
        V=V_next
        V_next=np.zeros((N_k))
    for ik, k in enumerate(k_grid): #ik is index number for the rows running from 0-->k_grid, while k in k_grid is an actual value
        for jk, kk in enumerate(k_grid): #jk is the index number for the columns
                Chi[ik,jk]=m[ik,jk]+beta*V[jk]
        V_next[ik] = np.nanmax(Chi[ik,:]) #value
        k_policy[ik] = np.argmax(Chi[ik,:])   #position: Policy Function capital
        c_policy[ik] =  y_func(k,1) +(1-delta)*k - k_policy[ik]   #Policy Functions consumption
      
    if np.allclose(V_next,V): #Break/Exit condition
        V_result=V_next
        break 
    else:
        i+=1
        
        
end_time = time.time() #Comparison of times: but are the same
T=qe.toc()    
print("total time taken this loop: ", end_time - start_time)   
print("total time taken this loop: ", T)  
print('total number of iterations: ', i)     
        
#Plotting Value Function iteration
fig, ax = plt.subplots()
ax.plot(k_grid,V_result, label='V(k)')
plt.title ('Value Function Iteration')
plt.legend()
plt.xlabel('Time') 
plt.show()

### Plotting Policy Functions
fig, ax = plt.subplots()
ax.plot(k_grid,c_policy, label='c(k)')
ax.plot(k_grid,k_policy, label='k(k)')
plt.title ('Policy Functions')
plt.legend()
plt.xlabel('Time') 
plt.show()

    
# =============================================================================
# 1.1 Value Function - bruce force iteration by QUANTECON Package- !running!
# =============================================================================      
V_guess=np.zeros((N_k))

def bellman_operator(V,return_policies=False):
    V_next = np.zeros((N_k))
    m = np.zeros((N_k,N_k)) #why define again?
    
    for ik, k in enumerate(k_grid):
        for igh, gh in enumerate(k_grid):
            ### Fill up return matrix: All possible choices of k' (k_grid) and labor (gh) per each state (k)
            m[ik,:] =  u(y_func(k,igh) +(1-delta)*k - k_grid) +beta*V
           
        V_next[ik] = np.nanmax(m[ik,:])     #wert   
        k_policy[ik] = np.argmax(m[ik,:])   #position
        c_policy[ik] =  y_func(k,1) +(1-delta)*k - k_policy[ik]         
    if return_policies==True:
        return V_next, k_policy, c_policy, m
    else:
        return V_next
    
#Compute fixed point in bellmann equation
qe.tic()
start_time = time.time()
V = qe.compute_fixed_point(bellman_operator, V_guess, max_iter=1000, error_tol=0.001, print_skip=20)
V, g_k, g_c, m = bellman_operator(V, return_policies=True)
qe.toc()
end_time=time.time()

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
i = 1 
j=0
V = np.zeros((N_k))
V_guess=np.zeros ((N_k))
V_next = np.zeros((N_k))
m = np.zeros((N_k,N_k))

while i<10000: 
    start_time = time.time()
    if i<=1: 
        V=V_guess 
    else: 
        V=V_next
        V_next=np.zeros((N_k))
    for ik, k in enumerate(k_grid):#ik is index number for the rows running from 0-->k_grid, while k in k_grid is an actual value
        
        for jk, kk in enumerate(k_grid): #jk is the index number for the columns
            if ik<jk: 
                m[ik,jk]= 0
            if kk<y_func(k,h)+(1-delta)*k:
                m[ik,jk] = u(y_func(k,h) +(1-delta)*k - kk)+beta*V[jk]
            else:
                m[ik,jk]= negativ #checking consumption >0
        V_next[ik] = np.nanmax(m[ik,:]) #value
        #Policy Function
        k_policy[ik] = np.argmax(m[ik,:])   #position
        c_policy[ik] =  y_func(k,1) +(1-delta)*k - k_policy[ik]   
      
    if np.allclose(V_next,V): 
        break 
    else:
        i+=1
        
end_time = time.time()     
print("total time taken this loop: ", end_time - start_time) 
print('total number of iterations: ', i)

fig, ax = plt.subplots()
ax.plot(k_grid,V_next, label='V(k) - monotonicity')
ax.plot(k_grid,V_result, label='V(k) - bruce force', ls='--')
plt.title ('Value Function Iteration')
plt.legend()
plt.xlabel('Time')  
plt.show()


# =============================================================================
# 1.3 Value Function - Concavity 
# =============================================================================
i = 1 
V = np.zeros((N_k))
V_guess=np.zeros ((N_k))
V_next = np.zeros((N_k))

#Value Iteration with Concavity
while i<10000: 
    start_time = time.time()
    if i<=1: 
        V=V_guess 
    else: 
        V=V_next
        V_next=np.zeros((N_k))
    for ik, k in enumerate(k_grid): #ik is index number for the rows running from 0-->k_grid, while k in k_grid is an actual value
        for jk, kk in enumerate(k_grid): #jk is the index number for the columns 
                if jk<2:
                    Chi[ik,jk]= m[ik,jk]+beta*V[jk]
                else:                                   #Concavity
                    jjk=jk-1
                    jjjk=jk-2
                    if Chi[ik,jjjk]>Chi[ik,jjk]:        #If the previous k_j-2 was greater than k_j-1 we can discard all following values
                        np.pad(Chi, (0,), 'constant', constant_values=0) #This fills up the array with zeros. A better option would evaluate directly the Array.
                    else:
                        Chi[ik,jk]= m[ik,jk]+beta*V[jk]
        V_next[ik] = np.nanmax(Chi[ik,:])            #Value Function
        k_policy[ik] = np.argmax(Chi[ik,:])          #Policy Function -position
        c_policy[ik] =  y_func(k,1) +(1-delta)*k - k_policy[ik]        #Policy Function 
      
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
ax.plot(k_grid,V_next, label='V(k) - bruce force', ls='--')
plt.show()

### Ploting Policy Functions
fig, ax = plt.subplots()
ax.plot(k_grid,c_policy, label='c(k)')
ax.plot(k_grid,k_policy, label='k(k)')
ax.legend()
plt.show()


# =============================================================================
# 1.5 Value Function - Concavity & Monotonicity
# =============================================================================
i = 1 
j=0
V = np.zeros((N_k))
V_guess=np.zeros ((N_k))
V_next = np.zeros((N_k))
m = np.zeros((N_k,N_k))

for ik, k in enumerate(k_grid): 
    for jk, kk in enumerate(k_grid): 
        if kk<y_func(k,h)+(1-delta)*k:
            m[ik,jk] = u(y_func(k,h) +(1-delta)*k - kk)    #checking consumption >0
        else:
            m[ik,jk] = negativ


while i<10000: 
    start_time = time.time()
    #(V_next-V_guess)>0.05:
#while i<1000:
    if i<=1: 
        V=V_guess 
    else: 
        V=V_next
        V_next=np.zeros((N_k))
    for ik, k in enumerate(k_grid): 
        for jk, kk in enumerate(k_grid): 
            if ik<jk:                               #Monotonicity
                Chi[ik,jk]= 0
            else: 
                if jk<2:
                    Chi[ik,jk]= m[ik,jk]+beta*V[jk]
                else:                                   #Concavity
                    jjk=jk-1
                    jjjk=jk-2
                    if Chi[ik,jjjk]>Chi[ik,jjk]:
                        np.pad(Chi, (0,), 'constant', constant_values=0)
                        #m.append()
                        #m.fill(0)
                    else:
                        Chi[ik,jk]= m[ik,jk]+beta*V[jk]
        V_next[ik] = np.nanmax(Chi[ik,:])            #Value Function
        k_policy[ik] = np.argmax(Chi[ik,:])          #Policy Function -position
        c_policy[ik] =  y_func(k,1) +(1-delta)*k - k_policy[ik]     #Policy Function 
      
    if np.allclose(V_next,V): 
        V_result=V_next
        break
    else:
        i+=1
               
end_time = time.time()     

print("total time taken this loop: ", end_time - start_time)   
print('total number of iterations: ', i)     
        
#Plotting graphs 
fig, ax = plt.subplots()
ax.plot(k_grid,V_result, label='V(k)')
ax.plot(k_grid,V_next, label='V(k) - bruce force', ls='--')
plt.show()

### Plots
fig, ax = plt.subplots()
ax.plot(k_grid,c_policy, label='c(k)')
ax.plot(k_grid,k_policy, label='k(k)')
ax.legend()
plt.show()





