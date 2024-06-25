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