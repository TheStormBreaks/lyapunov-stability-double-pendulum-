from __future__ import print_function
from scipy.integrate import odeint

import time
import math
import numpy as np
import pylab as py

#we import matplotlib.pyplot as plt.
from matplotlib import animation, rc
from IPython.display import HTML
from matplotlib import pyplot as plt


#mass of bob for pendulum 1 and pendulum 2 (in kg).
m1 = 5
m2 = 2

#length of string for pendulum 1 and oendulum 2 (in meters).
l1 = 2
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

#Conversion from polar coordinates to cartesian coordinates
#first pendulum
x1 = l1 * np.sin(u0)
y1 = -l1 * np.cos(u0)
#second pendulum
x2 = x1 + l2 * np.sin(u2)
y2 = y1 - l1 * np.cos(u2)


#MATPLOTLIB WITH PYPLOT

py.close('all')

#dot diagram

py.figure(1)
py.plot(x1, x2, '.', colour = '#0077BE', label = 'mass 1')
py.plot(x2, y2, '.', colour = '#f66338', label = 'mass 2')
py.legend()
py.xlabel('x (meters)')
py.ylabel('y (meters)')

#line graph 
py.figure(2)
py.plot(t,x2)
py.plot(t,y2)
py.legend()
py.xlabel('distance (meters)')
py.ylabel('time (seconds)')

fig = plt.figure()
ax = plt.axis(xlim = (-l1 - l1 - 0.5, l1 + l2 + 0.5), ylim = (-2.5, 1.5))

#line for pendulum1's path
line1, = ax.plot([], [], 'o-', colour = '#d2eeff', markersize = 12, markerfacecolor='#0077BE', 
                 lw = 2, markevery = 10000, markeredgecolour = 'k')
#line for pendulum2's path
line2, = ax.plot([], [], 'o-', colour = '#ffebd8', markersize = 12, markerfacecolor='#f66338', 
                 lw = 2, markevery = 10000, markeredgecolour = 'k')

#line3 and line4 are used to create lines in the animation plot
line3, = ax.plot([], [], color='k', linestyle='-', linewidth=2)
line4, = ax.plot([], [], color='k', linestyle='-', linewidth=2)

#This is a marker (specifically a circle marker 'o') used 
#to denote a specific point or event in the animation. 
line5, = ax.plot([], [], 'o', color='k', markersize=10)

#used to display the current time in the animation plot. 
time_template = 'Time = %.1f s'

#used to dynamically update and display the formatted time 
time_string = ax.text(0.05, 0.9, '', transform=ax.transAxes)

#The init() function serves as the initialization function for the animation.
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    line5.set_data([], [])
    time_string.set_text('')

    return line3, line4, line5, line1, line2, time_string

# sequence wise called animation function
def animate(i):
    #size of the motion trail. 
    # length of motion trail of weight 1
    trail1 = 6             
    # length of motion trail of weight 2
    trail2 = 8      

    # time step  
    dt = t[2]-t[1]  

    # marker + line of first weight
    line1.set_data(x1[i : max(1, i - trail1) : -1], y1[i : max(1, i - trail1) : -1])
    
    # marker + line of second weight  
    line2.set_data(x2[i : max(1, i - trail2) : -1], y2[i : max(1, i - trail2) : -1])

    #line connecting weight 2 to weight 1
    line3.set_data([x1[i], x2[i]], [y1[i], y2[i]])

    # line connecting the origin to the weight
    line4.set_data([x1[i], 0], [y1[i], 0])

    line5.set_data([0, 0], [0, 0])
    time_string.set_text(time_template % (i * dt))

    return line1, line2, line3, line4, line5, time_string

animate = animation.FuncAnimation(fig, animate, init_func = init, frames = Nt, 
                                  interval = 1000 * (t[2] - t[1]) * 0.8, blit = True)


#TO SAVE THE ANIMATION, EVERY RUN
#anim.save('double_pendulum_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
#anim.save('double_pendulum_animation.gif', fps=1.0/(t[2]-t[1]), writer='imagemagick')

plt.show()
