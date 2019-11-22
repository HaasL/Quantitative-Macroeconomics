#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 11:52:44 2019

@author: lillianhaas
"""
import numpy as np
import sympy as sy
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy.polynomial.polynomial as polinomial 
import numpy.polynomial.chebyshev as cheby
from mpl_toolkits.mplot3d import axes3d, Axes3D 
plt.style.use("ggplot")

# =============================================================================
# 1.2 Approximation RAMP Functions
# =============================================================================
'''
1. Define Grid 
2. Define Function 
3. Create Taylor Approximation  
3. Evaluation point along Taylor Approximation of the Ramp Function 
4. Plotting the graph'''

#1. Define the Grid 
x_lims = [-2, 6]
x1 = np.linspace(x_lims[0],x_lims[1],100)
y1 = []
x0=2 #point of evaluation

#2. Define the function 
# 1. Step Ramp Function
x=sy.Symbol('x', real=True)
def f(x):
    return 0.5*(x + abs(x))

#3. Create Taylor approximation 
def factorial(n): 
    if n<=0:
        return 1
    else: 
        return n*factorial(n-1)

def taylor(function,x0,n):
    i = 0 #i-th derivative of f evaluated at the point x
    p = 0
    while i <= n:
        p =p+ (function.diff(x,i).subs(x,x0))/(factorial(i)*(x-x0)**i)
        i += 1
    return p

#4. Approximation for degree of order [1,2,5,20]
x_lims = [-2, 6]
x1 = np.linspace(x_lims[0],x_lims[1],100)
y1 = []
y2 = []
y5 = [] 
y20 = []
approx = taylor(f(x),x0,1)
for k in x1:
    y1.append(approx.subs(x,k))
approx = taylor(f(x),x0,2)
for k in x1:
    y2.append(approx.subs(x,k))    
approx = taylor(f(x),x0,5)
for k in x1:
    y5.append(approx.subs(x,k))
approx = taylor(f(x),x0,20)
for k in x1:
    y20.append(approx.subs(x,k))
#5. Plotting the graph    
plt.plot(x1, f(x1), label='$ (x+|x|)/2 $')
plt.plot(x1,y1,label='order 1')
plt.plot(x1,y2,label='order 2')   
plt.plot(x1,y5,label='order 5')    
plt.plot(x1,y20,label='order 20')    
plt.xlim(x_lims)
plt.ylim([-2,6])
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.title('Taylor series approximation')
plt.show()

# =============================================================================
# 1.3.1 Approximation three Functions
# =============================================================================
'''
1. Define Grid 
2. Define all three functions
3. Define evenly spaced approximation nodes  
3. Interpolating cubic polynomials
4. Plotting the graph'''

#1. Define Grid
x_lims=[-1,1]
x_lims=[-1,1]
x1= np.linspace (x_lims[0],x_lims[1],11)
xs= np.linspace (x_lims[0],x_lims[1],101)

#Define all three functions
#First Function - 1/(1+25*x^2)
#Define the true function and the interpolation function
def h(x): 
    return 1/(1+25*x**2)
yh=h(x1)

#Second Function e^(1/x)
def j(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.exp(1/x)
yj=j(x1)
 
#Third Function (x+|x|)/)
def p(x):
    return 0.5*(x+abs(x))
yp=p(x1)

# Approximate up until 10 starting from 1 and using steps of 2
''' Try to make a shorter version...
order = [3, 5, 10]
ih= []
ihg= []
ij= []
ijg= []
ip= []
ipg= []

for k in order:
    ih.append(polinomial.polyfit(x1,yh,k))
    ihg.append(polinomial.polyval(xs,ih[-1]))
    ij.append(polinomial.polyfit(x1,yj,k))
    ijg.append(polinomial.polyval(xs,ij[-1]))
    ip.append(polinomial.polyfit(x1,yj,k))
    ipg.append(polinomial.polyval(xs,ij[-1]))'''

#Interpolating cubic polynomial, monomials of order 5 and 10    
interpolh=polinomial.polyfit(x1,yh,3)
interpol5h=polinomial.polyfit(x1,yh,5)
interpol10h=polinomial.polyfit(x1,yh,10)
interpolgh=polinomial.polyval(xs,interpolh)
interpol5gh=polinomial.polyval(xs,interpol5h)
interpol10gh=polinomial.polyval(xs,interpol10h)

#Interpolating cubic polynomial, monomials of order 5 and 10
interpolj=polinomial.polyfit(x1,yj,3)
interpol5j=polinomial.polyfit(x1,yj,5)
interpol10j=polinomial.polyfit(x1,yj,10)
interpolgj=polinomial.polyval(xs,interpolj)
interpol5gj=polinomial.polyval(xs,interpol5j)
interpol10gj=polinomial.polyval(xs,interpol10j)

#Interpolating cubic polynomial, monomials of order 5 and 10
interpolp=polinomial.polyfit(x1,yp,3)
interpol5p=polinomial.polyfit(x1,yp,5)
interpol10p=polinomial.polyfit(x1,yp,10)
interpolgp=polinomial.polyval(xs,interpolp)
interpol5gp=polinomial.polyval(xs,interpol5p)
interpol10gp=polinomial.polyval(xs,interpol10p)

#Figure plot Function - 1/(1+25*x^2)
plt.plot(x1, yh,'o', label= "nodes")
plt.plot(xs,h(xs), label= "True")
plt.plot (xs,interpolgh, label= "Cubic")
plt.plot (xs,interpol5gh, label= "Monomial5")
plt.plot (xs,interpol10gh, label= "Monomial10")
plt.legend()
plt.show()

#Figure plot Function - e^(1/x)
plt.plot(x1, yj,'o', label= "nodes")
plt.plot(xs,j(xs), label= "True")
plt.plot (xs,interpolgj, label= "Cubic")
plt.plot (xs,interpol5gj, label= "Monomial5")
plt.plot (xs,interpol10gj, label= "Monomial10")
plt.legend()
plt.show()

#Figure plot Function - (x+|x|)/2)
plt.plot(x1, yp,'o', label= "nodes")
plt.plot(xs,p(xs), label= "True")
plt.plot (xs,interpolgp, label= "Cubic")
plt.plot (xs,interpol5gp, label= "Monomial5")
plt.plot (xs,interpol10gp, label= "Monomial10")
plt.legend()
plt.show()

#Illustrating the error term 
plt.plot(xs,h(xs)-interpolgh,label='Cubic error') 
plt.plot(xs,h(xs)-interpol5gh,label='Monomial5 error')
plt.plot(xs,h(xs)-interpol10gh,label='Monomial10 error')
plt.ylim(-0.5,0.75)
plt.title('Approximation Errors - Runge Function')
plt.legend()
plt.show()

plt.plot(xs,j(xs)-interpolgj,label='Cubic error') 
plt.plot(xs,j(xs)-interpol5gj,label='Monomial5 error')
plt.plot(xs,j(xs)-interpol10gj,label='Monomial10 error')
plt.ylim(0,1)
plt.title('Approximation Errors - $(e^(1/x))$')
plt.legend()
plt.show()

plt.plot(xs,p(xs)-interpolgp,label='Cubic error') 
plt.plot(xs,p(xs)-interpol5gp,label='Monomial5 error')
plt.plot(xs,p(xs)-interpol10gp,label='Monomial10 error')
plt.ylim(-0.25,0.5)
plt.title('Approximation Errors - $((x+|x|)/2))$')
plt.legend()
plt.show()

# =============================================================================
# 1.3.2 Approximation three Functions
# =============================================================================
#Chebychev interpolation nodes
N=11
x1=np.cos((2*np.arange(1,N+1)-1)/(2*N)*np.pi)
x1=sorted (x1)
x1=np.asarray(x1)

#Interpolating Chebychev - 1/(1+25*x^2)
interpolh=polinomial.polyfit(x1,yh,3)
interpol5h=polinomial.polyfit(x1,yh,5)
interpol10h=polinomial.polyfit(x1,yh,10)
interpolgh=polinomial.polyval(xs,interpolh)
interpol5gh=polinomial.polyval(xs,interpol5h)
interpol10gh=polinomial.polyval(xs,interpol10h)

#Interpolating Chebychev - e^(1/x)
interpolj=polinomial.polyfit(x1,yj,3)
interpol5j=polinomial.polyfit(x1,yj,5)
interpol10j=polinomial.polyfit(x1,yj,10)
interpolgj=polinomial.polyval(xs,interpolj)
interpol5gj=polinomial.polyval(xs,interpol5j)
interpol10gj=polinomial.polyval(xs,interpol10j)

#Interpolating Chebychev - (x+|x|)/2)
interpolp=polinomial.polyfit(x1,yp,3)
interpol5p=polinomial.polyfit(x1,yp,5)
interpol10p=polinomial.polyfit(x1,yp,10)
interpolgp=polinomial.polyval(xs,interpolp)
interpol5gp=polinomial.polyval(xs,interpol5p)
interpol10gp=polinomial.polyval(xs,interpol10p)

#Plotting the figure
#Figure plot Function - 1/(1+25*x^2)
plt.plot(x1, yh,'o', label= "nodes")
plt.plot(xs,h(xs), label= "True")
plt.plot (xs,interpolgh, label= "Cubic")
plt.plot (xs,interpol5gh, label= "Monomial5")
plt.plot (xs,interpol10gh, label= "Monomial10")
plt.title('Chebychev Interpolation - Ramp function')
plt.legend()
plt.show()

#Figure plot Function - e^(1/x)
plt.plot(x1, yj,'o', label= "nodes")
plt.plot(xs,j(xs), label= "True")
plt.plot (xs,interpolgj, label= "Cubic")
plt.plot (xs,interpol5gj, label= "Monomial5")
plt.plot (xs,interpol10gj, label= "Monomial10")
plt.title('Chebychev Interpolation - $(e^(1/x))$')
plt.legend()
plt.show()

#Figure plot Function - (x+|x|)/2)
plt.plot(x1, yp,'o', label= "nodes")
plt.plot(xs,p(xs), label= "True")
plt.plot (xs,interpolgp, label= "Cubic")
plt.plot (xs,interpol5gp, label= "Monomial5")
plt.plot (xs,interpol10gp, label= "Monomial10")
plt.title('Chebychev Interpolation - $(x+|x|)/2)$')
plt.legend()
plt.show()

#Illustrating the error term 
plt.plot(xs,h(xs)-interpolgh,label='Cubic error') 
plt.plot(xs,h(xs)-interpol5gh,label='Monomial5 error')
plt.plot(xs,h(xs)-interpol10gh,label='Monomial10 error')
plt.ylim(-0.5,0.75)
plt.title('Approximation Errors by Cheb - Runge Function')
plt.legend()
plt.show()

plt.plot(xs,j(xs)-interpolgj,label='Cubic error') 
plt.plot(xs,j(xs)-interpol5gj,label='Monomial5 error')
plt.plot(xs,j(xs)-interpol10gj,label='Monomial10 error')
plt.ylim(0,1)
plt.title('Approximation Errors by Cheb  - $(e^(1/x))$')
plt.legend()
plt.show()

plt.plot(xs,p(xs)-interpolgp,label='Cubic error') 
plt.plot(xs,p(xs)-interpol5gp,label='Monomial5 error')
plt.plot(xs,p(xs)-interpol10gp,label='Monomial10 error')
plt.ylim(-0.25,0.5)
plt.title('Approximation Errors by Cheb  - $((x+|x|)/2))$')
plt.legend()
plt.show()

#################################
#Part 2 Chebychev polynomials
#################################

#Interpolation Function - 1/(1+25*x^2)
chebypolh=cheby.chebfit(x1,yh,3)
chebypol5h=cheby.chebfit(x1,yh,5)
chebypol10h=cheby.chebfit(x1,yh,10)
chebypolgh=cheby.chebval(xs,chebypolh)
chebypol5gh=cheby.chebval(xs,chebypol5h)
chebypol10gh=cheby.chebval(xs,chebypol10h)

#Interpolation  - e^(1/x)
chebypolj=cheby.chebfit(x1,yj,3)
chebypol5j=cheby.chebfit(x1,yj,5)
chebypol10j=cheby.chebfit(x1,yj,10)
chebypolgj=cheby.chebval(xs,chebypolj)
chebypol5gj=cheby.chebval(xs,chebypol5j)
chebypol10gj=cheby.chebval(xs,chebypol10j)

#Interpolation - (x+|x|)/2
chebypolp=cheby.chebfit(x1,yp,3)
chebypol5p=cheby.chebfit(x1,yp,5)
chebypol10p=cheby.chebfit(x1,yp,10)
chebypolgp=cheby.chebval(xs,chebypolp)
chebypol5gp=cheby.chebval(xs,chebypol5p)
chebypol10gp=cheby.chebval(xs,chebypol10p)

#Figure plot Function - 1/(1+25*x^2)
plt.plot(x1, yh,'o', label= "nodes")
plt.plot(xs,h(xs), label= "True")
plt.plot (xs,chebypolgh, label= "Cubic")
plt.plot (xs,chebypol5gh, label= "Monomial5")
plt.plot (xs,chebypol10gh, label= "Monomial10")
plt.title('Approximation by Chebychev polynomials - Runge function')
plt.legend()
plt.show()

#Figure plot Function - e^(1/x)
plt.plot(x1, yj,'o', label= "nodes")
plt.plot(xs,j(xs), label= "True")
plt.plot (xs,chebypolgj, label= "Cubic")
plt.plot (xs,chebypol5gj, label= "Monomial5")
plt.plot (xs,chebypol10gj, label= "Monomial10")
plt.title('Approximation by Chebychev polynomials - $(e^(1/x))$')
plt.legend()
plt.show()

#Figure plot Function - (x+|x|)/2
plt.plot(x1, yp,'o', label= "nodes")
plt.plot(xs,p(xs), label= "True")
plt.plot (xs,chebypolgp, label= "Cubic")
plt.plot (xs,chebypol5gp, label= "Monomial5")
plt.plot (xs,chebypol10gp, label= "Monomial10")
plt.title('Approximation by Chebychev polynomials - $((x+|x|)/2))$')
plt.legend()
plt.show()

#Illustrating the error term 
plt.plot(xs,p(xs)-chebypolgh,label='Cubic error-R') 
plt.plot(xs,p(xs)-chebypol5gh,label='Monomial5 error-R')
plt.plot(xs,p(xs)-chebypol10gh,label='Monomial10 error-R')
plt.ylim(-1,1)
#plt.title('Approximation Errors by Cheb  - Runge function')
#plt.legend()
plt.plot(xs,p(xs)-chebypolgj,label='Cubic error-e') 
plt.plot(xs,p(xs)-chebypol5gj,label='Monomial5 error-e')
plt.plot(xs,p(xs)-chebypol10gj,label='Monomial10 error-e')
#plt.ylim(-1,1)
#plt.ylim(-0.25,0.5)
#plt.title('Approximation Errors by Cheb  - $(e^(1/x))$ ')
#plt.legend()
plt.plot(xs,p(xs)-chebypolgp,label='Cubic error-Ramp') 
plt.plot(xs,p(xs)-chebypol5gp,label='Monomial5 error-Ramp')
plt.plot(xs,p(xs)-chebypol10gp,label='Monomial10 error-Ramp')
#plt.ylim(-1,1)
#plt.ylim(-0.25,0.5)
plt.title('Approximation Errors by Cheby - $((x+|x|)/2))$')
plt.title('Approximation Errors by Cheby polyinomials')
plt.legend()
plt.show()

# =============================================================================
# 2. Excersice 
# =============================================================================
# CES
#1. We have to create the Chebyshev nodes 
#2. Adopt these nodes to the 3-dimensional space 
#3. Evaluate our function at Chebyshev nodes + Chebyshev basis
#4. Compute Chebyshev coefficients along the Chebyshev basis

#Define the Function 
def f(alpha,sigma,k,h):
	return np.power((1-alpha)*k**((sigma-1)/sigma) + alpha*h**((sigma-1)/sigma),1/sigma)
alpha = .5
sigma = .25
a=0
b=10

# 1. Creating N=20 Chebychev nodes 
n = 20
m1=np.cos((2*np.arange(1,n+1)-1)/(2*n)*np.pi) 

#2. Shifting domain for k and h 
k1 =(b-a)/2*(1+m1)+a
h1 = k1
h1 = np.asarray(h1)
k1 = np.asarray(k1)
m1= np.asarray(m1)

#3. Evaluating true value of nodes at the function
nodes = f(alpha,sigma,k1,h1)
nodes = np.matrix(nodes)

#3. Setting the basis along the evaluated nodes
degree= 3
def T(d,x):
    r = []
    r.append(np.ones(len(x))) #1. polynomial along slides
    r.append(x) #2. polynomial along slides
    for i in range(1,d): 
        p = 2*x*r[i]-r[i-1]
        r.append(p)
    basis = np.matrix(r[d])
    return basis

#4. Calculating Chebychev coefficient along the formula presendet in slides 48
    #Here the computed nodes from above are weighted with the according polynomial defined in 3.
def coeff(b,nodes,d):
    theta=np.empty((d+1)*(d+1))
    theta.shape = (d+1,d+1)
    for i in range(d+1):
        for j in range(d+1):
            theta[i,j] = (np.sum(np.array(nodes)*np.array((np.dot(T(i,b).T,T(j,b)))))/np.array((T(i,b)*T(i,b).T)*(T(j,b)*T(j,b).T)))
    return theta

#5.The Approximation Function weighs along the chebyshev polynomials & coefficient each point of the graph 
def approximation(x,y,theta,d):
    f = []
    zh = (2*(y-a)/(b-a)-1)
    zk = (2*(x-a)/(b-a)-1)
    for u in range(d):
        for z in range(d):
                f.append(np.array(theta[u,z])*np.array((np.dot(T(u,zk).T,T(z,zh)))))
    approxsum = sum(f)
    return approxsum

# use m1 analogous to z_k in the picture: 
theta = coeff(m1, nodes, degree)
true = f(alpha, sigma, k1[:, None], h1[None, :])
approx = approximation(k1, h1, theta, degree)
error = f(alpha, sigma, k1[:, None], h1[None, :]) - approximation(k1, h1, theta, degree)

#Plotting the figure 
X, Y = np.meshgrid(k1,h1)
fig = plt.figure()
#ax = Axes3D(fig) 
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, true)
ax.set(title='True')
ax.set_xlabel('K')
ax.set_ylabel('H')
ax.set_zlabel('f(K,H)')
ax.set_xlim(0,10)
ax.set_ylim(0,10)
ax.set_zlim(0,10)
plt.show()

###Try with degree 6 
degree=6
theta = coeff(m1, nodes, degree)
true = f(alpha, sigma, k1[:, None], h1[None, :])
approx = approximation(k1, h1, theta, degree)
error = f(alpha, sigma, k1[:, None], h1[None, :]) - approximation(k1, h1, theta, degree)

#Plotting the figure 
X, Y = np.meshgrid(k1,h1)
fig = plt.figure()
#ax = Axes3D(fig) 
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, true)
ax.set(title='True')
ax.set_xlabel('K')
ax.set_ylabel('H')
ax.set_zlabel('f(K,H)')
ax.set_xlim(0,10)
ax.set_ylim(0,10)
ax.set_zlim(0,10)
plt.show()


    