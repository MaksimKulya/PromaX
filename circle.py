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


def circle(Nx,Ny,X,Y,r):
    
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)

    circ = np.zeros(shape=(Nx, Ny))

    for i in range(Nx):
        for j in range(Ny):
            if x[i]**2+y[j]**2<r**2:
                circ[i,j]=1
            else:
                circ[i,j]=0
                
    return circ








