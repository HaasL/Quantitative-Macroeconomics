#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 12:31:16 2019

@author: lillian
"""

import numpy as np 
import matplotlib.pyplot as plt 
#import collections, numpy

# =============================================================================
# Excercise 1.1 -Generating the data set along Multivariate Normal Distribution 
# =============================================================================
N=10000000
k_i= np.array((N))
z_i=np.array((N))

covar1= [[1, 0], [0, 1]]
covar2= [[1, 0.5], [0.5, 1]]
covar3= [[1, -0.5], [-0.5, 1]]
i=1

Runs=[covar1, covar2, covar3]

for cov in Runs:

#Generating random values of normal variat function
    mean = [1, 1]
    cov = [[1, 0], [0, 1]]
    k, z = np.random.multivariate_normal(mean, cov, 10000000).T

#Plotting the graph in logs & levels
    k_i=k[:100000]
    z_i=z[:100000] #to many point to plot 
    plt.plot(k_i, z_i,'')
    plt.axis('equal')
    plt.title('Joint distribution ln(k) & ln(z)') #plot in logs, what are levels?
    plt.ylabel('ln(z)')
    plt.xlabel('ln(k)')
    plt.show()

# =============================================================================
# Excercise 1.2 - Computing output 
# =============================================================================
#saving generated distribution values 
    k.T
    z.T
    k_i=np.exp(k)
    z_i=np.exp(z)
    gamma= 0.6
    
    y_i=k_i**gamma*z_i

# =============================================================================
# Excercise 1.3 - Finding optimal output
# =============================================================================
    z_s=np.zeros(10000000)
    k_o=np.zeros(10000000)

    z_s=(z_i[0]/z_i)**(1/(gamma-1))    #computing first z-elements
    a=sum(k_i)/(sum(z_s))#**1/(gamma-1))   #computing first k[0]

    k_o=a*z_s
    k_o[0]=a    
    print ('The sum of optimal capital is',sum(k_o), 'while the random pick was is',sum(k_i)) 
# =============================================================================
# Excercise 1.4 - Compare "data" to equilibria 
# =============================================================================

    plt.plot(k_o,ls='--')
    plt.plot(k_i,color='c')
    plt.ylabel('')
    plt.xlabel('')
    plt.show()
    
    print ('Mean k_i',np.mean(k_i), 'in Run',i)
    print ('Mean k_o',np.mean(k_o),'in Run',i)

# =============================================================================
# Excercise 1.5 - Outout gains
# =============================================================================
    y_o=k_o**gamma*z_i
    Y_o=sum(y_o)
    Y_i=sum(y_i)

    gains=1

#Gains from computation 
    gains=((Y_o/Y_i)-1)*100 # first result: -100
  
    print('Gains from optimal capital implementation for covariance',i,'achieves output gains of', gains)
    i+=1
# =============================================================================
# Excercise 1.6 - Re-do with new Correlations
# =============================================================================
#See the loop above 
    
# =============================================================================
# Excercise 2 - Re-do with new gamma
# =============================================================================
#New 
gamma= 0.8      #NEW GAMMA
covar1= [[1, 0], [0, 1]]
covar2= [[1, 0.5], [0.5, 1]]
covar3= [[1, -0.5], [-0.5, 1]]
i=1

Runs=[covar1, covar2, covar3]


for cov in Runs:

#Generating random values of normal variat function
    mean = [1, 1]
    k2, z2 = np.random.multivariate_normal(mean, cov, 10000000).T


# Excercise 2.1.2 - Computing output 
#saving generated distribution values 
    k2.T
    z2.T
    k_i=np.exp(k2)
    z_i=np.exp(z2)        
    y_i=k_i**gamma*z_i
    

# Excercise 2.1.3 - Finding optimal output
    z_s=np.zeros(10000000)
    k_o=np.zeros(10000000)

    z_s=(z_i[0]/z_i)**(1/(gamma-1))    #computing first z-elements
    a=sum(k_i)/(sum(z_s))#**1/(gamma-1))   #computing first k[0]

    k_o=a*z_s
    k_o[0]=a    
    print ('The sum of optimal capital is',sum(k_o), 'while the random pick was is',sum(k_i)) 

# Excercise 2.1.4 - Compare "data" to equilibria 
    
    print ('Mean k_i',np.mean(k_i), 'in Run',i)
    print ('Mean k_o',np.mean(k_o),'in Run',i)

# Excercise 2.1.5 - Outout gains
    y_o=k_o**gamma*z_i
    Y_o=sum(y_o)
    Y_i=sum(y_i)

    gains=1

#Gains from computation 
    gains=(Y_o/Y_i-1)*100 # first result: -100
  
    print('Gains from optimal capital implementation for covariance',i,'achieves output gains of', gains)
    i+=1

# =============================================================================
# Excercise 3.1 - Correlation, Variance of 10.000 random observations
# =============================================================================
#Pay attention here we consider the ln of k and ln z 
#while i<=10000: 
k_c=np.random.choice(k,10000)
z_c=np.random.choice(z,10000)

print ('Variance of 10.000 k is', np.var(k_c))
print ('Variance of 10.000 z is', np.var(z_c))
print ('Co-variance of 10.000 k is', np.cov(k_c))


# =============================================================================
# Excercise 3.2 - Redo 1.1-1.5
# =============================================================================
gamma= 0.6
k_c.T
z_c.T
k_1=np.zeros(10000)
del y_i 
del k_i 
del z_i 
del y_o 

k_i=np.exp(k_c)
z_i=np.exp(z_c)
#Output based on random k
y_i=k_i**gamma*z_i    

# Excercise 2.1.3 - Finding optimal capital trajectory 
z_s=np.zeros(10000)
k_o=np.zeros(10000)     #introducing optimal capital

z_s=(z_i[0]/z_i)**(1/(gamma-1))    #computing first z-elements
a=sum(k_i)/(sum(z_s))#**1/(gamma-1))   #computing first k[0]

k_o=a*z_s
k_o[0]=a    
print ('The sum of optimal capital is',sum(k_o), 'while the random pick was is',sum(k_i)) 

# Excercise 2.1.4 - Compare "data" to equilibria 
plt.plot(k_o,ls='--')
plt.plot(k_i,color='c')
plt.ylabel('k')
plt.title('Comparison of trajectory paths of capital for sample size 10.000')
plt.show()
    
print ('Mean k_i',np.mean(k_i))
print ('Mean k_o',np.mean(k_o))

# Excercise 2.1.5 - Outout gains
y_o=k_o**gamma*z_i 
Y_o=sum(y_o)
Y_i=sum(y_i)

#Gains from computation 
gains=(Y_o/Y_i-1)*100
print ('Gains of 10.000 sample selection',gains)
# =============================================================================
# Excercise 3.3 - Repeat 1000 times and show histogram, median
# =============================================================================
i=0
gains=[]

while i<1000: 
    k_1=np.zeros(10000)
    k_c=np.random.choice(k,10000)
    z_c=np.random.choice(z,10000)

    k_c.T
    z_c.T

    k_i=np.exp(k_c)
    z_i=np.exp(z_c)
    #Output based on random k
    k_1=k_i**gamma
    y_i=np.multiply(z_i,k_1)
    
    # Excercise 2.1.3 - Finding optimal capital trajectory 
    z_s=np.zeros(10000)
    k_o=np.zeros(10000)     #introducing optimal capital
    z_s=(z_i[0]/z_i)**(1/(gamma-1))    #computing first z-elements
    a=sum(k_i)/(sum(z_s))#**1/(gamma-1))   #computing first k[0]
    k_o=a*z_s
    k_o[0]=a    

    # Excercise 2.1.5 - Outout gains
    k_1=k_o**gamma
    y_o=np.multiply(z_i,k_1)

    #Gains from computation 
    gains.append(((sum(y_o)/sum(y_i))-1)*100)
    i+=1    

plt.hist(gains)
plt.title('Distribution of output gains for 1.000 random samples of size 10.000')
plt.ylabel('k')
np.median(gains)

# =============================================================================
# Excercise 3.4 - Probability to get close of 10%
# =============================================================================
#Assuming 138 is optimum
#For 10% divergence
low=138/100*95
high=138/100*105
j=0
k=0
i=0
b=0

while i<1000:
    if gains[i]>=low:
        b+=1 #values greater equal low
        if gains[i]<=high:
            j+=1 #values in boundries
    else: 
        k+=1 #outside boundries
    i+=1    
print ('Probability to be close to th optimum is', j/1000)    

#For 10% divergence in each direction
low=138/100*90
high=138/100*110          
j=0
k=0
i=0
b=0
while i<1000:
    if gains[i]>=low:
        b+=1 #values greater equal low
        if gains[i]<=high:
            j+=1 #values in boundries
    else: 
        k+=1 #outside boundries
    i+=1    
print ('Probability to be close to th optimum is', j/1000)   

# =============================================================================
# Excercise 3.5 - Repeat for different sample sizes 
# =============================================================================
gamma=0.6
ratio1= 100
ratio2= 1000
ratio3= 100000


Runs=[ratio1, ratio2, ratio3]

for N in Runs:
    i=0
    gains=[]
    
    while i<1000: 
        k_1=np.zeros(N)
        k_c=np.random.choice(k,N)
        z_c=np.random.choice(z,N)
    
        k_c.T
        z_c.T
        k_i=np.exp(k_c)
        z_i=np.exp(z_c)
        #Output based on random k
        y_i=k_i**gamma*z_i
        
        # Excercise 2.1.3 - Finding optimal capital trajectory 
        z_s=np.zeros(N)
        k_o=np.zeros(N)     #introducing optimal capital
        z_s=(z_i[0]/z_i)**(1/(gamma-1))    #computing first z-elements
        a=sum(k_i)/(sum(z_s))#**1/(gamma-1))   #computing first k[0]
        k_o=a*z_s
        k_o[0]=a    
        
        # Excercise 2.1.5 - Outout gains
        y_o=k_o**gamma*z_i
        #Gains from computation 
        gains.append((sum(y_o)/sum(y_i)-1)*100)
        i+=1    
        ja=0+1/1000
    print('Result for sample size',N)    
    plt.hist(gains)
    plt.title('Result for sample size')
    plt.ylabel('k')
    plt.legend(gains)
    print('The mean for sample size',N,' output gains is',np.median(gains))
    