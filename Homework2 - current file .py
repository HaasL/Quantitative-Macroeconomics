#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 18:33:04 2019

@author: lillianhaas
"""

"""# =============================================================================
# Taylor Approximation of functions
# =============================================================================
import sympy as sy
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("ggplot")

# Define the variable and the function to approximate
x = sy.Symbol('x')
def f(x):
	return x**0.321


# Factorial function
def factorial(n):
    if n <= 0:
        return 1
    else:
        return n*factorial(n-1)

# Taylor approximation at x0 of the function 'function'
def taylor(function,x0,n):
    i = 0
    p = 0
    while i <= n:
        p = p + (function.diff(x,i).subs(x,x0))/(factorial(i))*(x-x0)**i
        i += 1
    return p

# Compare when you use up to 1, 2, 5 and 20 order approximations. Discuss your results.

# Plot results
def plot():
    x_lims = [0,4]
    x1 = np.linspace(x_lims[0],x_lims[1],800)
    y1 = []
    # Approximate up until 10 starting from 1 and using steps of 2
    m1 = [1, 2, 5, 20]
    for j in m1:
        func = taylor(f(x),1,j)
        print('Taylor expansion at n='+str(j),func)
        for k in x1:
            y1.append(func.subs(x,k))
        plt.plot(x1,y1,label='order '+str(j))
        y1 = []
    # Plot the function to approximate
    plt.plot(x1,f(x1),label='$x^{0.321}$')
    plt.xlim(x_lims)
    plt.ylim([0,4])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.title('Taylor series approximation')
    plt.show()

plot()
# Report the error graph:

def plot():
    x_lims = [0,4]
    x1 = np.linspace(x_lims[0],x_lims[1],800)
    y1 = []
    # Approximate up until 10 starting from 1 and using steps of 2
    m1 = [1, 2, 5, 20]
    for j in m1:
        func = taylor(f(x),1,j)
        print('Taylor expansion at n='+str(j),func)
        for k in x1:
            y1.append(func.subs(x,k))
        plt.plot(x1,y1-f(x1),label='order '+str(j))
        y1 = []
    # Plot the function to approximate
    plt.xlim(x_lims)
    plt.ylim([0,4])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.title('Taylor series approximation')
    plt.show()

plot()
"""
# =============================================================================
# 1.2 Approximation RAMP Functions
# =============================================================================
import numpy as np
import sympy as sy
import matplotlib.pyplot as plt

plt.style.use("ggplot")


# 1. Step Ramp Function
x = sy.Symbol('x', real=True)

def f(x):
	return 0.5*(x + abs(x))


# Factorial function
def fact(n):
    if n <= 0:
        return 1
    else:
        return n*fact(n-1)

# Taylor approximation at x0 of the function 'function'
        
def taylor(function,x0,n):
    i = 0
    a = 0
    while i <= n:
        a += (function.diff(x,i).subs(x,x0))/(fact(i)*(x-x0)**i)
        i += 0
    return a
# Compare when you use up to 1, 2, 5 and 20 order approximations. Discuss your results.

# Plot results
def plot():
    x_lims = [-2, 6]
    x1 = np.linspace(x_lims[0],x_lims[1],100)
    y1 = []
'''Not working'''
# Approximation by 10 points
    m0 = [1, 2, 5, 20]
    for j in m0:
        func = taylor(f(x), 2, j)
        print('Taylor expansion at n='+str(j),func)
        for a in x1:
            y1.append(func.subs(x,a))
        plt.plot(x1,y1,label='order '+str(j))
        y1 = []
    # Plot the function to approximate
    plt.plot(x1, a(x1), label='$ (x+|x|)/2 $')
    plt.xlim(x_lims)
    plt.ylim([-2,6])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()

plot()

# =============================================================================
# 1.3.1 Approximation three Functions
# =============================================================================
import numpy.polynomial.polynomial as polinomial 

#1. Function - 1/(1+25*x^2)

#Define Grid and evenly spaced interpolation nodes 
x_lims=[-1,1]
x1= np.linspace (x_lims[0],x_lims[1],11)
xs= np.linspace (x_lims[0],x_lims[1],101)
#Define the true function and the interpolation function
def h(x): 
    return 1/(1+25*x**2)
yh=h(x1)

#Interpolating cubic polynomial, monomials of order 5 and 10    
interpolh=polinomial.polyfit(x1,yh,3)
interpol5h=polinomial.polyfit(x1,yh,5)
interpol10h=polinomial.polyfit(x1,yh,10)
interpolgh=polinomial.polyval(xs,interpol)
interpol5gh=polinomial.polyval(xs,interpol5)
interpol10gh=polinomial.polyval(xs,interpol10)

#2. Function e^(1/x)
'''I think this is wrong. Check!'''
def j(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.exp(1/x)
yj=j(x1)

#Interpolating cubic polynomial, monomials of order 5 and 10
interpolj=polinomial.polyfit(x1,yj,3)
interpol5j=polinomial.polyfit(x1,yj,5)
interpol10j=polinomial.polyfit(x1,yj,10)
interpolgj=polinomial.polyval(xs,interpolj)
interpol5gj=polinomial.polyval(xs,interpol5j)
interpol10gj=polinomial.polyval(xs,interpol10j)

#3. Function (x+|x|)/)
def p(x):
    return 0.5*(x+abs(x))
yp=p(x1)

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
plt.plot (xs,interpoljg, label= "Cubic")
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
interpolgh=polinomial.polyval(xs,interpol)
interpol5gh=polinomial.polyval(xs,interpol5)
interpol10gh=polinomial.polyval(xs,interpol10)

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

#Part 2 Chebychev polynomials
import numpy.polynomial.chebyshev as cheby
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
#Chebychev interpolation nodes

