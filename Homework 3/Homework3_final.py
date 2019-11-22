#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 16:55:53 2019

@author: lillian
"""
import numpy as np 
import matplotlib.pyplot as plt
import sympy as sy
from scipy.optimize import fsolve
import pandas as pd 

# =============================================================================
# Excercise 1
# =============================================================================
#Given data by PS3 
theta=0.67
ht=0.31
yss1=1.0        	
      	
# =============================================================================
# 1.1 & 1.2 Analytically solved for steady state in t=1 for z and z_new=2*2
# =============================================================================
#By computation as shown in pdf document:
z=1.62968            	
beta= 0.98039  
delta=0.0625  
k1=4.0 
k2=8.0     	

# =============================================================================
# 1.3 Transition path from steady state 1 to steady state 2
# =============================================================================
#1. Situation of steady state
#ct=0.955 #First guess for c after shock
ct=1.0
#ct=1.02   #interesting to see how sensitive the curve reacts to a slight change in c
c=[]
k=[]
y1=[]
inv=[]
kt=4
k.append(kt)
c.append(ct)
inv.append(delta*kt)
y1.append(yss1)

#2. Shock occurs
z*=2

#3. The market reacts. We observe a transition of the market
#I had a iteration loop for the initial guess on ct but it was not working

while (k2-kt>0.1):
    ktn=kt**(1-theta)*(z*ht)**(theta)+(1-delta)*kt-ct #ktn is k_{t+1}
    ctn=beta*ct*((1-theta)*(z*ht)**(theta)*ktn**(-1*theta)+(1-delta)) #Euler Equation
    inv.append(ktn-(1-delta)*kt)
    ytn=ktn**(1-theta)*(z*ht)**theta
    y1.append(ytn)
    k.append(ktn)
    c.append(ctn)
    kt=ktn.real
    ct=ctn.real

#Plotting the figure  
plt.figure()
plt.subplot(221)
plt.plot(c)
plt.title('Consumption')
plt.xlabel('Time')
plt.subplot(222)
plt.plot(k)
plt.title('Capital')
plt.xlabel('Time')
plt.subplot(223)
plt.plot(y1)
plt.title('Output')
plt.xlabel('Time')
plt.subplot(224)
plt.plot(inv)
plt.title('Investment')
plt.xlabel('Time')
plt.subplots_adjust(top=2, bottom=0.08, left=0, right=2, hspace=0.3, wspace=0.2)
plt.show()    

#Plotting the figure
b=len(k)
print ("It takes",b,"time periods to converge to the new steady state")

print (k,c,inv,y1)

# =============================================================================
# 1.4 Transition path from steady state 1 to steady state 2 and back to steady state 1
# =============================================================================
'''1. Considering transition path as in 1.3 up to t=10 
2. Guessing the reacting of the agents behaviour by consumption in t=11
3. Computing the transition to steady state 
4. Checking if the transition is reached for the initial guess of consumption 
5. Adjusting the guess of consumption in t=11 in order to approach the first steady state
'''
#Shock in t=10
z/=2

#1. Step Same transition as seen in 1.3
cs=[]
ks=[]
ys=[]
invs=[]

cs.extend(c[:9])
ks.extend(k[:10])
ys.extend(y1[:10])
invs.extend(inv[:10])
#First guess for c after shock
#ct=c[-1]-.1
ct=0.905
cs.append(ct)
kt=ks[-1]

#Transition path
#I had a iteration loop for the initial guess on ct but it was not working
while (kt-k1 > 0.001):
    ktn=kt**(1-theta)*(z*ht)**(theta)+(1-delta)*kt-ct #ktn is k_{t+1}
    ctn=beta*ct*((1-theta)*(z*ht)**(theta)*ktn**(-1*theta)+(1-delta)) #Euler Equation
    invs.append(ktn-(1-delta)*kt)
    ytn=ktn**(1-theta)*(z*ht)**theta
    ys.append(ytn)
    ks.append(ktn)
    cs.append(ctn)
    kt=ktn.real
    ct=ctn.real

#4. Plotting Transition path 
plt.figure()
plt.subplot(221)
plt.plot(cs)
plt.title('Consumption')
plt.xlabel('Time')
plt.subplot(222)
plt.plot(ks)
plt.title('Capital')
plt.xlabel('Time')
plt.subplot(223)
plt.plot(ys)
plt.title('Output')
plt.xlabel('Time')
plt.subplot(224)
plt.plot(invs)
plt.title('Investment')
plt.xlabel('Time')
plt.subplots_adjust(top=2, bottom=0.08, left=0, right=2, hspace=0.3, wspace=0.2)
plt.show() 

b=len(ks)
print ("It takes",b,"times periods to converge to the new steady state")

print (ks,cs,invs,ys)

# =============================================================================
# Excercise 2
# =============================================================================

# =============================================================================
# 2.1 General equilibrium in closed economies
# =============================================================================
#Comment by TA:Regarding the last question on optimal taxation, it will depend on the welfare function of the social planner.
#Since the exercise doesn’t specify any form, you have freedom to choose and discuss.
'''
1. 6 Unknown parameters are: w,r,h and c respectively for (h and low)
2.Solving differential equations analytically for Langragian for latter 4 parameters 
and the profit function for w and r. 
3. By setting a system of equations Python is able to solve it.
4. The system of equations shows for each closed economy (A,B) the optimal solution 
'''
#Assuming capital shares for low and high qualified workers
kl_a=1.0
kh_a=1.0
kl_b=1.0
kh_b=1.0
#Given parameters (switched for low and high skill!)
kappa=5.0
nu=1.0
sigma=0.8
prodl_a=0.5 
prodl_b=2.5
prodh_a=5.5 
prodh_b=3.5
z=1.0
theta=0.6 #the value of lambda should approx clean the government budget constraint so you do not have to worry for that. 
k_2=kl_a+kh_a #equals 2
lambda_a=0.95
lambda_b=0.84
phi=0.2

'''Notation: 
    r= a1=0
    w=a2 =1
    hl=a3=2
    hh=a4=3
    cl=a5=4
    ch=a6=5'''

#Country A
#Note: the end _a indicates country A, the end _b indicates country B

def a(x): 
    a1=(-x[0]+(1-theta)*z*(x[2]*prodl_a+x[3]*prodh_a)**(theta)*k_2**(-theta))
    a2=(-x[1]+(theta)*z*k_2**(1-theta)*(x[2]*prodl_a+x[3]*prodh_a)**(theta-1))
    a3=(-x[2]+((1-phi)*lambda_a*(1/kappa)*x[4]**(-1*sigma)*(x[1]*prodl_a)**(1-phi))**(nu/(1+nu*phi)))
    a4=(-x[3]+((1-phi)*lambda_a*(1/kappa)*x[5]**(-1*sigma)*(x[1]*prodh_a)**(1-phi))**(nu/(1+nu*phi)))
    a5=(-x[4]+lambda_a*(x[1]*x[2]*prodl_a)**(1-phi)+x[0]*kl_a**prodl_a)
    a6=(-x[5]+lambda_a*(x[1]*x[3]*prodh_a)**(1-phi)+x[0]*kh_a**prodh_a)
    return (a1, a2, a3, a4, a5, a6)
solution1a= fsolve(a,[1,1,1,1,1,1])
print (solution1a)

print ('Country A, interest rate r:',+solution1a[0])
print ('Country A, wages:',+solution1a[1])
print ('Country A, labor - low skill:',+solution1a[2])
print ('Country A, labor - high skill:',+solution1a[3])
print ('Country A, consumption - low skill:',+solution1a[4])
print ('Country A, consumption - high skill:',+solution1a[5])


#Country B

def b(x): 
    b1=-x[0]+((1-theta)*z*(x[2]*prodl_b+x[3]*prodh_b)**(theta)*k_2**(-theta))
    b2=-x[1]+((theta)*z*k_2**(1-theta)*(x[2]*prodl_b+x[3]*prodh_b)**(theta-1))
    b3=-x[2]+(((1-phi)*lambda_b*(1/kappa)*x[4]**(-1*sigma)*(x[1]*prodl_b)**(1-phi))**(nu/(1+nu*phi)))
    b4=-x[3]+(((1-phi)*lambda_b*(1/kappa)*x[5]**(-1*sigma)*(x[1]*prodh_b)**(1-phi))**(nu/(1+nu*phi)))
    b5=-x[4]+(lambda_b*(x[1]*x[2]*prodl_b)**(1-phi)+x[0]*kl_b**prodl_b)
    b6=-x[5]+(lambda_b*(x[1]*x[3]*prodh_b)**(1-phi)+x[0]*kh_b**prodh_b)
    return (b1,b2,b3,b4,b5,b6)
solution1b= fsolve(b,[1,1,1,1,1,1])
print (solution1b)    

print ('Country B, interest rate r:',+solution1b[0])
print ('Country B, wages:',+solution1b[1])
print ('Country B, labor - low skill:',+solution1b[2])
print ('Country B, labor - high skill:',+solution1b[3])
print ('Country B, consumption - low skill:',+solution1b[4])
print ('Country B, consumption - high skill:',+solution1b[5])
    
    
# =============================================================================
# 2.1 General equilibrium in closed economies
# =============================================================================
'''
1. 16 Unknown parameters are: w,r,h,c and k all for country A and B. 
Parameters h,c,k exist also for each type (high and low) individually 
2.Differential equations for w and r are solved from the Production function. 
The other paramters are solved for by the utlity function and its budget contraint for consumption 
3. Python solves for the system of equations 
4. The solutions shows the optimum values in the GE equilibrium for these two countries A & B
'''
#Note: the end _a indicates country A, the end _b indicates country B
#kl_a=1.0
#kh_a=1.0
#kl_b=1.0 
#kh_b=1.0 
k_2=4 
#16 Unknowns to find the equilibrium 
'''
#Notation
    r_a= a1=0
    r_b=a2=1
    w_a=a3=2
    w_b=a4=3
    hl_a=a5=4
    hh_a=a6=5
    hl_b=a7=6
    hh_b=a8=7 
    cl_a=a9=8
    ch_a=a10=9
    cl_b=a11=10
    ch_b=a12=11
    kl_a=a13=12
    kh_a=a14=13
    kl_b=a15=14
    kh_b=a16=15'''

#System of equations to solve for the second optimization problem

def f(x): 
    f1=-x[0]+(1-theta)*z*(prodh_a*x[5]+prodl_a*x[4])**(theta)*k_2**(-1*theta) 
    f2=-x[1]+(1-theta)*z*(prodl_bx[6]+prodh_b*x[7])**(theta)*k_2**(-1*theta)
    f3=-x[2]+(theta)*z*k_2**(1-theta)*(x[5]+x[4])**(theta-1)
    f4=-x[3]+(theta)*z*k_2**(1-theta)*(x[6]+x[7])**(theta-1)
    f5=-x[4]+((1-phi)*lambda_a*(1/kappa)*x[8]**(-1*sigma)*(x[2]*prodl_a)**(1-phi))**(nu/(1+nu*phi))
    f6=-x[5]+((1-phi)*lambda_a*(1/kappa)*x[9]**(-1*sigma)*(x[2]*prodh_a)**(1-phi))**(nu/(1+nu*phi))
    f7=-x[6]+((1-phi)*lambda_b*(1/kappa)*x[10]**(-1*sigma)*(x[3]*prodl_b)**(1-phi))**(nu/(1+nu*phi))
    f8=-x[7]+((1-phi)*lambda_b*(1/kappa)*x[11]**(-1*sigma)*(x[3]*prodh_b)**(1-phi))**(nu/(1+nu*phi))
    f9=-x[8]+lambda_a*(x[2]*x[4]*prodl_a)**(1-phi)+x[0]*kl_a**prodl_a+x[1]*(k_2-kl_a)
    f10=-x[9]+lambda_a*(x[2]*x[5]*prodh_a)**(1-phi)+x[0]*kh_a**prodh_a+x[1]*(k_2-kh_a)
    f11=-x[10]+lambda_b*(x[3]*x[6]*prodl_b)**(1-phi)+x[1]*kl_b**prodl_b+x[0]*(k_2-kl_b)
    f12=-x[11]+lambda_b*(x[3]*x[7]*prodh_b)**(1-phi)+x[1]*kh_b**prodh_b+x[0]*(k_2-kh_b)
    f13=-x[12]+(x[1]/(x[0]*prodl_a))**(1/(prodl_a-1))
    f14=-x[13]+(x[1]/(x[0]*prodh_a))**(1/(prodh_a-1))
    f15=-x[14]+(x[0]/(x[1]*prodl_b))**(1/(prodl_b-1))
    f16=-x[15]+(x[0]/(x[1]*prodh_b))**(1/(prodh_b-1))
    return (f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16)

solution2=fsolve(f,[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])

print ('Country A+B, interest rate r in country A:',+solution2[0])
print ('Country A+B, interest rate r in country B:',+solution2[1])
print ('Country A+B, wages w in country A:',+solution2[2])
print ('Country A+B, wages rate w in country B:',+solution2[3])
print ('Country A+B, labor - low skill in country A:',+solution2[4])
print ('Country A+B, labor - high skill in country A:',+solution2[5])
print ('Country A+B, labor - low skill in country B:',+solution2[6])
print ('Country A+B, labor - high skill  in country B:',+solution2[7])
print ('Country A+B, consumption - low skill in country A:',+solution2[8])
print ('Country A+B, consumption - high skill in country A:',+solution2[9])
print ('Country A+B, consumption - low skill in country A:',+solution2[10])
print ('Country A+B, consumption - high skill in country B:',+solution2[11])
print ('Country A+B, capital - low skill in country A:',+solution2[12])
print ('Country A+B, capital - high skill in country A:',+solution2[13])
print ('Country A+B, capital - low skill in country B:',+solution2[14])
print ('Country A+B, capital - high skill in country B:',+solution2[15])
    
