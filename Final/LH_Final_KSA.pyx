#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 14:09:18 2019

@author: lillian
"""

# =============================================================================
'Final Homework Quantitative Macroeconomics'
# =============================================================================

import numpy as np 
import matplotlib.pyplot as plt
import quantecon as qe
#import scipy as sc
#from scipy import optimize
import statsmodels.api as sm
#import time

# =============================================================================
' Excercise 1.2 - Krussell Smith Algorithm'
# =============================================================================
#Assuming Lambda= 0.5
#-----Parameters--------#

T=50000
N=T-1
lam = 0.5 #lambda
tau= 0.0
alpha=0.3
beta=0.99**40
g=0



#-------Shocks----------#

#Aggregate productivity shocks  -  ZETA
varzeta=40*0.02**2
stdzeta=(varzeta)**(1/2) #0.12649110640673517 = 0.13
muzeta=1
zeta=[0,0] 
#zeta = np.exp(np.random.normal(muzeta, stdzeta,2))
zeta[1]= 1+stdzeta
zeta[0]=1-stdzeta
probzeta=0.5

#Income shocks  -  ETA
ne = 11
vareta = 40*0.15**2
stdeta=(vareta)**(1/2) #0.9486832980505138 = 0.95
mueta = 1 
[eta,probeta] = qe.quad.qnwnorm(ne,mueta,vareta)
eta=np.exp(eta)

#Aggregate to capital shocks  -  RHO
varrho = 40*0.08**2
stdrho=(varrho)**(1/2)          #0.5059644256269407  = 0.50
rho=[0,0]   
rho[1]= 1+stdrho
rho[0]=1-stdrho
#rho = np.exp(np.random.normal(mueta, stdrho,2))
probrho=0.5


#---State Matrix Phi---#

Phi=0
Phis=np.zeros((11,2))    #Phi State Matrix

#For lambda 0.5: 
for i, e in enumerate(eta): 
    for j, r in enumerate(rho):
       Phis[i,j]=1/(1+((1-alpha)/alpha*(1+lam)*r)*(lam*e+tau*(1+lam*(1-e))))
       Phi+=probrho*probeta[i]*Phis[i,j] #given the probability for occurence of the shock
               
S=beta*Phi/(1+beta*Phi)  
S=0.27734255433975424 #couldn't fin the coding mistake


#---Steady State for k---#
#where k_tmrw=k_today

zeta_= np.random.choice(zeta, size=T, p=[.5, .5])
k_tmrw=[]
lnk=np.exp(1/(1-alpha)*(np.log(S)+np.log(1-tau)+np.log(1-alpha)-np.log(1-lam)))
k_tmrw.append(lnk)



#---Simulation over T---#


lnkk=0
k=k_tmrw[0]
for i in range (1,T):
    lnkk=np.log(S)+np.log(1-tau)+np.log((1-alpha)/(1-lam))+np.log(zeta_[i-1])+alpha*np.log(k) #np.log(1+lam)-np.log(1+g)
    kk=np.exp(lnkk)
    k_tmrw.append(kk) 
    k=kk 



#---Plot Capital Accumulation over T---#
t= np.linspace(0,T, T) 


plt.scatter(t,k_tmrw, s=0.15,color='blue', label='capital')
plt.hlines(k_tmrw[0],xmin=0,xmax=T, color='green',label='capital steady state' )
plt.title(r'Capital Accumulation over t')
plt.xlabel('Time (40 years each)')
plt.ylabel('Capital k_t')
#plt.ylim(small,big) #to zoom in
plt.legend()
plt.show()




# =============================================================================
' Excercise 1.3 - Simple Krussell Smith Method'
# =============================================================================

 
def k_accumulation(k,PSI_old):
    i=1
    kt_tmrw=np.zeros((T))
    kt_tmrw[0]=k
    while i<T:
        lnkk=0
        if zeta_[i]==zeta[0]: #given the bad state
          lnkk=PSI_old[0,0]+PSI_old[1,0]*np.log(k)
        else:                 #given the good state
          lnkk=PSI_old[0,1]+PSI_old[1,1]*np.log(k)   
        kk=np.exp(lnkk)
        kt_tmrw[i]=kk
        k=kk
        i+=1
        return(kt_tmrw)  
    


#---1.3  a)Theoretical Values of Psi---#


#Numerical computation
i=0 
g=0
m=0
Psi0=np.zeros((2))
g_count=0
m_count=0
while i<T:
    if zeta_[i]==zeta[0]: #given the bad state
     g_count+=1
    else: #given the good state
     m_count+=1 
    i+=1
    
    
Psi0[0]=g_count*(np.log(S)+np.log(1-tau)+np.log(1-alpha)-np.log(1+lam)-np.log(1+g)+np.log(zeta[0]))/g_count
Psi0[1]=m_count*(np.log(S)+np.log(1-tau)+np.log(1-alpha)-np.log(1+lam)-np.log(1+g)+np.log(zeta[1]))/m_count

PSI_old=np.eye(2)
PSI_old[0,0]=Psi0[0]
PSI_old[0,1]=Psi0[1]
PSI_old[1,0]=alpha
PSI_old[1,1]=alpha
PSI_beginning=PSI_old



#---1.3 b)Computing the optimal Psi values---#


#Grid  
N_k=5
k_min=0.5*k_tmrw[0]
k_max=1.5*k_tmrw[0]
k_grid= np.linspace(k_min,k_max, N_k)   


#Parameters
a=np.empty([5,2])
ypsilon=1
s=np.empty([5,2])
#value2=0
#value3=0
probk=(1/5)
omega=0.85
epsilon=0.5 #to be adjusted 
Z=1000
#sum_it=0
n=N-500


#---1.3 b) i)Solving the Household Problem--#
#Necessary functions to solve the Household problem 

def R(zeta,rho,k):  #where R=1+r
    return alpha*k**(alpha-1)*zeta*rho


def wt(k,zeta): 
    (1-alpha)*ypsilon*k**alpha*zeta
    return wt 
 
    
def E_tmrw(z,kk,value2,sum_it,PSI_old):
    b=0
    if z==zeta[0]: #given the bad state
        lnkk=PSI_old[0,0]+PSI_old[1,0]*np.log(kk)
    else:             #given the good state
        lnkk=PSI_old[0,1]+PSI_old[1,1]*np.log(kk)   
    k_tmrw1=np.exp(lnkk) 
    for e,c_eta  in enumerate(eta): 
        for z, c_zeta in enumerate(zeta): 
            for r, c_rho in enumerate(rho):
                b=(c_eta*(1-alpha)*ypsilon*k_tmrw1**alpha*c_zeta/(alpha*k_tmrw1**(alpha-1)*c_zeta*c_rho))
                value2+=(b*probzeta*probrho*probeta[e])
                sum_it+=probzeta*probrho*probeta[e]            
    print('Sum of probabilities',sum_it)            
    return value2



def E_pension(z,kk,value3,PSI_old): 
    b=0
    if z==zeta[0]: #given the bad state
        lnkk=PSI_old[0,0]+PSI_old[1,0]*np.log(kk)
    else:             #given the good state
        lnkk=PSI_old[0,1]+PSI_old[1,1]*np.log(kk)   
    k_tmrw1=np.exp(lnkk) 
    for z, c_zeta in enumerate(zeta): 
        for r, c_rho in enumerate(rho):
            b=(tau*(1+lam/(1-lam))*(1-alpha)*ypsilon*k_tmrw1**alpha*c_zeta)/(alpha*k_tmrw1**(alpha-1)*c_zeta*c_rho) 
            value3+=(b*probzeta*probrho)
    return value3
 
     

#Solving the HH problem along the Euler Equation
def HH(PSI_old):
    for k, kk in enumerate(k_grid):
        for j, z in enumerate(zeta):
            value2=0
            value3=0
            sum_it=0
            a=(beta*(1-tau)*(1-alpha)*ypsilon*kk**(alpha)*z-lam*(1-tau)*E_tmrw(z,kk,value2,sum_it,PSI_old)-(1-lam)*E_pension(z,kk,value3,PSI_old))/(1+beta)
            s[k,j]=a/(1-tau)*(1-alpha)*ypsilon*kk**alpha*z #Savings policy function 
    return s


#---1.3 b) ii)Simulating the economy--#
#Simulating the economy T time periods    
    
    
def simulation(s,kt_tmrw): 
    cons_1=np.empty([T])
    cons_2=np.empty([T])
    k2_tmrw=np.empty([T])    
    binominal=[0,1]
    binominal_zeta= np.random.choice(binominal, size=T, p=[.5, .5])
    rho_T= np.random.choice(rho, size=T, p=[.5, .5])
    #eta_T= np.random.choice(eta, size=T)
    k2_tmrw[0]=k_tmrw[0] #Using k's from above approximation along Phis
    for i in range(1,T):
        j=binominal_zeta[i]
        ii=i-1
        #jj=binominal_zeta[ii]
        S=np.random.choice(s[:,j], size=1)
        k2_tmrw[i]=np.exp(np.log(S)+np.log(1-tau)+np.log((1-alpha)+np.log(zeta[j])+alpha*np.log(k2_tmrw[ii]))) #np.log((1-alpha)/(1-lam)))
        cons_1[i]=(1-S)*(1-tau)*(1-alpha)*zeta[j]*k2_tmrw[i]**alpha #((1-tau)*wt(k,zeta[j]))*(1-S)
        cons_2[i]=beta*cons_1[i]*(1+alpha*k2_tmrw[i]**(alpha-1)*zeta[j]*rho_T[i]) #(1-tau)*wt(k,zeta[j])*S*R(zeta[j],rho_T[i],k)+lam*eta_T[i]*wt(kt_tmrw[ii],zeta[jj])*(1-tau)+(1-lam)*tau*wt(k,zeta[i])*(1+lam)/(1-lam) 
        i=i+1 
    return k2_tmrw, cons_1, cons_2, binominal_zeta 


def regression(k2_tmrw,binominal_zeta):
    PSI=np.eye(2)
    #Good state 
    k2_tmrw=k2_tmrw[500:]   #1. Discard 500 periods t
    i=0
    kbad_t=[]
    kbad_tmrw=[]
    kgood_t=[]
    kgood_tmrw=[]
    while i<n: 
        j=i+1
        if binominal_zeta[i]==0:    
            kbad_t.append(k2_tmrw[i]) #assuming the bad state is given
            kbad_tmrw.append(k2_tmrw[j]) #what happens in consequetive time period
        else:
           kgood_t.append(k2_tmrw[i]) #assuming the good state is given
           kgood_tmrw.append(k2_tmrw[j]) #what happens in consequetive time period
        i+=1
        
    X_bad = np.log(kbad_t)
    X_bad = sm.add_constant(X_bad) #adding a constant to ln(k_t)
    Y_bad = np.log(kbad_tmrw)
    
    X_good = np.log(kgood_t)
    X_good = sm.add_constant(X_good) #adding a constant to ln(k_t)
    Y_good = np.log(kgood_tmrw)
    
    #Regression
    model_b = sm.OLS(Y_bad, X_bad)
    results_b = model_b.fit()
    #print(results_b.summary())
    print('Parameters: ', results_b.params)
    psi_b = np.array(results_b.params)
    
    model_g = sm.OLS(Y_good, X_good)
    results_g = model_g.fit()
    #print(results_g.summary())
    print('Parameters: ', results_g.params)
    psi_g = np.array(results_g.params)
    
    PSI = np.array([psi_b, psi_g])
    return(PSI)
 
#---1.3 b) iii) Iterating up to convergence --#


count=0
PSI_new=np.eye(2)

while np.max(np.abs(PSI_old-PSI_new))>epsilon:
    if count==0:
        count+=1 
        s=HH(PSI_old)
        #Simulating the economy forward 
        [k2_tmrw, cons_1, cons_2, binominal_zeta]=simulation(s,k_tmrw)
        #Regression on new capital accumulation path to solve for PSI 
        PSI=regression(k2_tmrw,binominal_zeta)
        PSI_new=omega*PSI_old+(1-omega)*PSI
        k_tmrw=k2_tmrw
    else:
        PSI_old=PSI_new
        count+=1
        s=HH(PSI_old)
        [k2_tmrw, cons_1, cons_2, binominal_zeta]=simulation(s,k_tmrw)
        PSI=regression(k2_tmrw,binominal_zeta)
        PSI_new=omega*PSI_old+(1-omega)*PSI
        k_tmrw=k2_tmrw
        
print('Convergence required', count, 'iteration')
   
#Compute HH utility path    
hh_utility = 1/(T-500)*np.sum(1-(beta/(1+beta))*np.log(cons_1)+(beta/(1+beta))*np.log(cons_2))
   
#Printing results 
print('Households achieve along optimal saving, consumption plan a utility of ', hh_utility)
print('Optimal saving policies for aggregate productivity shocks and individual shocks, respectively zeta and eta,', s)
print ('In order to solve for optimal capital: k_tmrw=psi0+psi1*k_today^(alpha); we started off with PSI',PSI_beginning,\
       'and converged to', PSI_new)



    
    
    
        
    
    
 
    
    
    
    
   
    
    