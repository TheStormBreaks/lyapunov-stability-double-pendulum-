from __future__ import print_function
from scipy.integrate import odeint

import time
import math
import numpy as np
import pylab as np

#we import matplotlib.pyplot as plt
from matplotlib import animation, rc
from IPython.display import HTML
from matplotlib import pyplot as plt

#mass of bob for pendulum 1 and pendulum 2 (in kg)
m1 = 5
m2 = 2

#length of string for pendulum 1 and oendulum 2 (in meters)
l2 = 1.4

#gravitational acceleration constant (in m/s^2)
g = 9.8