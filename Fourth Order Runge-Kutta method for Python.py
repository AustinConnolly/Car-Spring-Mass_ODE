 # Fourth Order Runge-Kutta method for Python
# CNNAUA001
# 23/04/2020

import numpy as np
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import Symbol


## Runge-Kutta Method for finding the stady states of Y(t) 

def RK4(v): # This is the Runge-Kutta method
    
    ## Initialised Variables: ----------------------------------------------------------------------
    
    a = 5*10**(-2)
    lamda = 10
    c = 2*10**3
    k = 7*10**4
    m = 800    
    omega = (2 * np.pi * v) / lamda
    omega_0 = np.sqrt(k / m) 
    gamma = c / m
    beta = ( (omega_0**2 - omega**2)**2 + (omega * gamma)**2 )      
    Y = np.zeros(601)
    h = 0.1
    z = 1 #-gamma/2 +1+((omega*(a * omega*(gamma**2 - omega_0**2)+a*omega_0**4 ))/beta)
    y = 1 #(1-(omega**3*gamma*a)/beta)
    t = 10/gamma
    i = 0
    Steadystate = np.zeros(500)
    maxa = np.zeros(500)
    ##----------------------------------------------------------------------------------------------
    
    ## Definig Equations: --------------------------------------------------------------------------
    
    def f(t,y,z): # This is the function z=dy/dt
        return z

    def g(t,y,z): # This is the function dz/dt = d^2y/dt^2
    
        a = 5*10**(-2)
    
        lamda = 10
    
        c = 2*10**3
    
        k = 7*10**4
    
        m = 800
    
        omega = (2 * np.pi * v) / lamda
    
        omega_0 = np.sqrt(k / m) 
    
        gamma = c / m    
    
        n = gamma*a*omega*np.cos(omega*t) + omega_0**2*a*np.sin(omega*t) - gamma*z - (omega_0**2)*(y)
        return n    
    ##----------------------------------------------------------------------------------------------
    
    ## Solving the Equation usind RK4 method: ------------------------------------------------------
    
    while t<60:
        
        ## Determining k1,k2,k3,k4:
    
        k1 = h * f(t,y,z)
        l1 = h * g(t,y,z)
        k2 = h * f(t + h/2, y + k1/2, z + l1/2)
        l2 = h * g(t + h/2, y + k1/2, z + l1/2)
        k3 = h * f(t + h/2, y + k2/2, z + l2/2)
        l3 = h * g(t + h/2, y + k2/2, z + l2/2)
        k4 = h * f(t+h, y+k3, z+l3)
        l4 = h * g(t+h, y+k3, z+l3)
        
        
        ## Plugging into equation:
        
        
        
        y = y + 1/6 * (k1 + 2*k2 + 2*k3 + k4) # Updates y-values
        
        Y[i] = y # Stores the y-values moved from below y=y...
        
        z = z + 1/6 * (l1 + 2*l2 + 2*l3 + l4) # Updates z-values
        
        i+=1 # Incriments the array Y
        
        t+=h # Increases the t-value
        
    ## Getting the Steady State from the function Y(t):
    for i in range(100,600):
        Steadystate[i-100] = Y[i]        
    maxa=np.amax(Steadystate)
    return maxa # Returns the Steady States of the function Y(t)
##--------------------------------------------------------------------------------------------------

'''
#---------------------------------------------------------------------------------------------------
## Initialising variables for finding the maximum amplitudes:
c=372
v = np. linspace(0.5,40,c)
Amax=np.zeros(c)
j=0
Max_v=np.zeros(c)
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
## Getting the maximum amplitudes

for i in v:
    Amax[j] = RK4(i)
    
    if Amax[j]>=np.amax(Amax):
        Max_v[j]=i  # Stores the velocities to later get the max velocity
    j+=1
print(np.amax(Max_v))
print(np.amax(Amax))
#print(np.amax(RK4(14.511022044088175)))

'''
#---------------------------------------------------------------------------------------------------
'''
## Question 2 A(v) vs v:

a = 5*10**(-2)

lamda = 10

c = 2*10**(3)

k = 7*10**4

m = 800

#v = 14.51102204488175

#v = Symbol( 'v' )

omega = (2 * np.pi * v) / lamda

omega_0 = np.sqrt(k / m) 

gamma = c / m

alpha = (np.sqrt(4 * omega_0**2 - gamma**2) ) / 2

beta = ( (omega_0**2 - omega**2)**2 + (omega * gamma)**2 )

A = -(omega**3 * gamma * a) / beta 

B = (a*omega**2*(gamma**2 - omega_0**2) + a* omega_0**4)/ beta

Atilde = np.sqrt( A**2 + B**2 )


#---------------------------------------------------------------------------------------------------
##Plottig the graphs:

plt.plot(v,Amax,'b-',label='Numerical')
plt.plot(v,Atilde,'r-',label='Analytical')
#plotother graph
xmin,xmax,ymin,ymax = plt.axis([0,40,0,0.2])

plt.tick_params(direction = 'in', top = 'on', right = 'on')
plt.xlabel('velocity [m/s]')
plt.ylabel(' max. amplitude [m]')
#plt.title('A graph showing the analytical and numerical amplitude graphs',size=10)
#plt.title('A graph showing the analytical maximum amplitude vs velocities of the car',size=10)
plt.title('A graph showing the numerical maximum amplitude vs velocities of the car',size=10)
plt.grid()#color='k',linestyle='-',linewidth=1*10**(-10)
plt.legend()
plt.show()

## Calculations:

s=(np.amax(Atilde))
print(s)
#print(solve((0.19514367204477992)**2 - A**2 - B**2,v  )) # solves the equation 0.19515564402464997)**2 - A**2 - B**2 = 0, for v


m=14.553908355795148 +gamma/2
mn=14.553908355795148-gamma/2

w=14.6076794657766 +gamma/2
wm=14.6076794657766-gamma/2

print(m,mn,w,wm)
'''
v=[15.8576794657766,15.803908355795148,13.3576794657766,13.303908355795148]

print(RK4(v[3]))
print('Done')#print(Atilde,np.amax(Amax),'done')
#---------------------------------------------------------------------------------------------------

print(100-(14.553908355795148/14.6076794657766)*100)
print((0.033/0.19515564402464997)*100)