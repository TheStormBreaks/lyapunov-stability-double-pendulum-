from __future__ import print_function
from scipy.integrate import odeint

import time
import math
import numpy as np
import pylab as np

#we import matplotlib.pyplot as plt.
from matplotlib import animation, rc
from IPython.display import HTML
from matplotlib import pyplot as plt


#mass of bob for pendulum 1 and pendulum 2 (in kg).
m1 = 5
m2 = 2

#length of string for pendulum 1 and oendulum 2 (in meters).
l2 = 1.4

#gravitational acceleration constant (in m/s^2).
g = 9.8

#initial conditions.
u0 = [-np.pi/2.2, 0, np.pi/1.8, 0]

# u[0] = angle of first pendulum
# u[1] = angular velocity of the first pendulum
# u[2] = angle of the second pendulum
# u[3] = angle valocity of the second pendulum

#Final time, Simulation time = 0 to tfinal.
tfinal = 25.0
Nt = 751
t = np.linspace (0, tfinal, Nt)
#the final time, in this case, sets the duration of the simulation
#for the double pendulum system
#It determines how long the simulation will run and for how many 
#time steps it will calculate the motion of the pendulum system. 
#In this programe, the final time is set to 25.0 seconds, meaning,
#the simulation will run for 25 seconds only. Then, it resets. 



#Differential equations describing the system 
def double_pendulum(u, t, m1, m2, l1, l2, g):
    #dU = derivatives
    # u = variables
    # p = parameters
    # t = time variable

    du = np.zeros(4)

    #intermediate variables
    c = np.cos(u[0] - u[1])
    s = np.sin(u[0] - u[2])

    du[0] = u[1]
    du[1] = (m2 * g * np.sin(u[2]) * c - m2 * s * (l1 * c * u[1] ** 2 + l2 * u[3] ** 2) - 
             (m1 + m2) * g * np.sin(u[0])) / (l1 * (m1 + m2 * s ** 2))
    du[2] = u[3]
    du[3] = ((m1 + m2) * (l1 * u[1] ** 2 * s - g * np.sin(u[2]) + g * np.sin(u[0]) * c) + 
             m2 * l1 * u[3] ** 2 * s * c) / (l2 * (m1 + m2 * s ** 2))
    
    return du

sol = odeint(double_pendulum, u0, t, args=(m1, m2, l1, l2, g))

# sol[:,0] = u1 = Θ_1
# sol[:,1] = u2 = ω_1
# sol[:,2] = u3 = Θ_2
# sol[:,3] = u4 = ω_2

#theta_1
u0 = sol[:, 0]

#omega_1
u1 = sol[:, 1]

#theta_2
u2 = sol[:, 2]

#omega_2
u3 = sol[:, 3]