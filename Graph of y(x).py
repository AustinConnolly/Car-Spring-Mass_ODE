# Plots graph of y(x) in OD Assignment 1

import matplotlib.pyplot as plt
import numpy as np


##Initialising Variables------------------------------------------------------------------------------------------------------------------------------

a = 5*10**(-2)

lamda = 10

c = 2*10**3

k = 7*10**4

m = 800

v = 12 #unknown currently

omega = (2 * np.pi * v) / lamda

omega_0 = np.sqrt(k / m) 

gamma = c / m

alpha = (np.sqrt(4 * omega_0**2 - gamma**2) ) / 2

beta = ( (omega_0**2 - omega**2)**2 - (omega * gamma)**2 )

x = np.linspace(-5.0,30.0)


## Graphing-------------------------------------------------------------------------------------------------------------------------------------------

y = a * np.sin((omega / v) * x)


## Plotting-------------------------------------------------------------------------------------------------------------------------------------------

plt.plot(x,y, 'b-') # plot time on x-axis, angle on y-axis  'b-' specifies a blue line
    
plt.xlabel('x')
    
plt.ylabel('y')
    
plt.title('A graph showing the Energy of the pendulum vs time')
    
#xmin, xmax, ymin, ymax = plt.axis ([0, 10, 30, 110])
    
plt.tick_params(direction = 'in', top = 'on', right = 'on')

plt.show()