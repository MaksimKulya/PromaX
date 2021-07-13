import matplotlib.pyplot as plt
import math
import cmath
import pylab
from matplotlib import mlab
from PIL import Image
import numpy as np
import scipy.constants
import numba as nb

def Gauss(Nx,Ny,X,Y,rho):

    x = np.linspace(-X/2, X/2-X/Nx, Nx)
    y = np.linspace(-Y/2, Y/2-Y/Ny, Ny)

    Gauss2d = np.zeros(shape=(Nx, Ny))

    for i in range(Nx):
        for j in range(Ny):
            Gauss2d[i, j] = np.exp(-2*(x[i]**2 + y[j]**2)/rho**2)

    return Gauss2d