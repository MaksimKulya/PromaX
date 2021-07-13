import matplotlib
import matplotlib.pyplot as plt
import math
import cmath
# import pylab
# from matplotlib import mlab
from PIL import Image
import numpy as np
import scipy.constants
import time



zz1 = np.linspace(0,50,11)
zz2 = np.linspace(50,500,10)
zz3 = np.linspace(500,2000,11)

zz = np.concatenate((zz1[0:-1], zz2[0:-1], zz3))

for i in zz:
    print(i)