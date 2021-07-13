from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import random
import matplotlib
import matplotlib.pyplot as plt
import math
import cmath
# import pylab
# from matplotlib import mlab
from PIL import Image
import numpy as np
import scipy.constants
import scipy
import time

fig = pyplot.figure()
ax = Axes3D(fig)

x=np.linspace(-5,5,100)
y=np.linspace(-5,5,100)
z=np.linspace(0,100,100)

ax.scatter(x, y, z)
pyplot.show()