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
def Gshift(GGin,nu,c,z,N0):
    
    GGout = np.zeros(shape=(GGin.shape[0],GGin.shape[1],GGin.shape[2]), dtype=nb.types.complex64)
    
    for k in range(GGin.shape[0]):
        for i in range(GGin.shape[1]):
            for j in range(GGin.shape[2]):
                GGout[k,i,j] = GGin[k,i,j]*np.exp(1j*2*pi*nu[k]*N0*z/c)
                
    return GGout


