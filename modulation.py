import matplotlib.pyplot as plt
import math
import cmath
import pylab
from matplotlib import mlab
from PIL import Image
import numpy as np
import scipy.constants
import numba as nb

pi=math.pi

@nb.njit
def modulation(Nx,Ny,am,ph,nu,n,c,h):

    AM = np.zeros(shape=(nu.shape[0], Nx, Ny))
    PH = np.zeros(shape=(nu.shape[0], Nx, Ny))
    G_object = np.zeros(shape=(nu.shape[0], Nx, Ny), dtype=nb.types.complex64)

    for k in range(nu.shape[0]):
        AM[k, :, :] = am
        PH[k, :, :] = ph * (2*pi*nu[k]*(n-1)*h/c)
        G_object[k, :, :] = AM[k, :, :] * np.exp(1j * PH[k, :, :])

    return G_object
