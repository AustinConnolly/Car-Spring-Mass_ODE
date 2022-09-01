#
#
#

import matplotlib.pyplot as plt
import numpy as np
from sympy.solvers import solve
from sympy import Symbol



##Initialising Variables------------------------------------------------------------------------------------------------------------------------------

#v = Symbol('v')

a = 5*10**(-2)

lamda = 10

c = 2*10**(3)

k = 7*10**4

m = 800

#v = np.linspace(-50,50,600)

v=[15.8576794657766,15.803908355795148,13.3576794657766,13.303908355795148]

omega = (2 * np.pi * v[3]) / lamda

omega_0 = np.sqrt(k / m) 

gamma = c / m

alpha = (np.sqrt(4 * omega_0**2 - gamma**2) ) / 2

beta = ( (omega_0**2 - omega**2)**2 + (omega * gamma)**2 )

A = -(omega**3 * gamma * a) / beta 

B = (a*omega**2*(gamma**2 - omega_0**2) + a* omega_0**4)/ beta

Atilde= np.sqrt( A**2 + B**2 )


print(Atilde)

#print(solve((0.19515564402464997)**2 - A**2 - B**2,v  )) # solves the equation (0.19515564402464997)**2 - A**2 - B**2 = 0, for v
#0.19515564402464997
#print(max(Atilde))
print('done')
##Plotting--------------------------------------------------------------------------------------------------------------------------------------------
'''
plt.plot(v,Atilde, 'b-') # plot time on x-axis, angle on y-axis  'b-' specifies a blue line
    
plt.xlabel('v [km/h]')
    
plt.ylabel('Amplitude [km]')
    
#plt.title('A graph showing the Energy of the pendulum vs time')
    
#xmin, xmax, ymin, ymax = plt.axis ([0,1,0,0.022])
    
plt.tick_params(direction = 'in', top = 'on', right = 'on')

plt.grid(True)



plt.show()'''