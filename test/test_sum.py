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


a1=np.ones(shape=(50,2,2))

a2=a1.sum(axis=(1,2))

Q=np.zeros(shape=(60))

numin=39;numax=99
for i in range(numin,numax):
    Q[i-numin]=1

