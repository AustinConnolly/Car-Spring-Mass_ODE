# Getting the maximum amplitudes
# graph Y(t) with different v's and record the maximum amplitude, from this make a table A(v) vs v to see if it matches the result we got analytically 


import matplotlib.pyplot as plt
import numpy as np


##Initialising Variables------------------------------------------------------------------------------------------------------------------------------

a = 5*10**(-2)

lamda = 10

c = 2*10**3

k = 7*10**4

m = 800



#omega = (2 * np.pi * v) / lamda

#omega_0 = np.sqrt(k / m) 

gamma = c / m

#alpha = (np.sqrt(4 * omega_0**2 - gamma**2) ) / 2

#beta = ( (omega_0**2 - omega**2)**2 + (omega * gamma)**2 )

#A = -(omega**3 * gamma * a) / beta 

#B = (a*omega**2*(gamma**2 - omega_0**2) + a* omega_0**4)/ beta

#N = 40  # number of iterations

Amax = np.zeros(400) # max amplitude of Y(t)

#t = np.linspace(10/gamma,50,600)

#Y=np.zeros(N)

#v=np.zeros(N)

#t=np.zeros(N)

#omega=np.zeros(N)

#beta=np.zeros(N)

## Functions -----------------------------------------------------------------------------------------------------------------------------------------

#Atilde = np.sqrt( A**2 + B**2 )

#Y = np.exp((-gamma/2)*t)*( np.cos(alpha * t)) + np.exp((-gamma/2)*t)*np.sin(alpha * t)- ((omega**3*gamma*a)/beta) * np.cos(omega*t)+ ((a*omega**2*(gamma**2-omega_0**2) + a * omega_0**4)/beta) * np.sin(omega*t) # Y(t) function

#print(max(Y))

#Q = np.linspace()
'''
for i in range(N):
    
    v[i] = i
    
    omega[i] = (2 * np.pi * v[i]) / lamda
    
    omega_0 = np.sqrt(k / m) 
    
    gamma = c / m
    
    alpha = (np.sqrt(4 * omega_0**2 - gamma**2) ) / 2
    
    beta[i] = ( (omega_0**2 - omega[i]**2)**2 + (omega[i] * gamma)**2 )
    
    #while t<1000:
    
    t = np.linspace(10/gamma,500)
    
    Y[i] = np.abs(np.amax(np.exp((-gamma/2)*t)*( np.cos(alpha * t)) + np.exp((-gamma/2)*t)*np.sin(alpha * t)- ((omega[i]**3*gamma*a)/beta[i]) * np.cos(omega[i]*t)+ ((a*omega[i]**2*(gamma**2-omega_0**2) + a * omega_0**4)/beta[i]) * np.sin(omega[i]*t))) # Y(t) function
    
    
       
    #Amax[i] = np.amax(Y)
    
    i+=1
 
#print(Amax)

#v=np.linspace(1,50,600)

plt.plot(v,Y,'b-')'''


t = np.linspace(10/gamma,40,400)

def A(v):
    omega = (2 * np.pi * v) / lamda
    
    omega_0 = np.sqrt(k / m) 
    
    gamma = c / m
    
    alpha = (np.sqrt(4 * omega_0**2 - gamma**2) ) / 2
    
    beta = ( (omega_0**2 - omega**2)**2 + (omega * gamma)**2 ) 
    
    Y = np.exp((-gamma/2)*t)*( np.cos(alpha * t)) + np.exp((-gamma/2)*t)*np.sin(alpha * t)- ((omega**3*gamma*a)/beta) * np.cos(omega*t)+ ((a*omega**2*(gamma**2-omega_0**2) + a * omega_0**4)/beta) * np.sin(omega*t)

    return np.amax(Y)

vel=np.zeros(400)
Q=np.linspace(0.5,40,400)
p=0
for i in Q:
    Amax[p]=A(i)
   
    if Amax[p]==np.amax(Amax):
        vel[p]=i
    i+=1
    p+=1     
print(np.amax(vel))
plt.plot(Q,Amax)
print(np.amax(Amax))
##Plotting--------------------------------------------------------------------------------------------------------------------------------------------

#plt.plot(t,Y, 'b-') # plot time on x-axis, angle on y-axis  'b-' specifies a blue line
    
plt.xlabel('v')
    
plt.ylabel('A(v)')
    
#plt.title('A graph showing the Energy of the pendulum vs time')
    
#xmin, xmax, ymin, ymax = plt.axis ([0, 100, 0, 0.210])
    
plt.tick_params(direction = 'in', top = 'on', right = 'on')

plt.show()




