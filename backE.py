import matplotlib.pyplot as plt
import math
import cmath
import pylab
from matplotlib import mlab
from PIL import Image
import numpy as np
import scipy.constants
import time
import Gshift

def backE(GG,N,nu,c,z,N0):

    GGshift=Gshift.Gshift(GG,nu,c,z,N0)    

    GG_ansig = np.pad(GGshift, ((N//2+1, N//2-GGshift.shape[0]-1), (0, 0), (0, 0)), mode='constant')
    
    EE = np.zeros(shape=(GG_ansig.shape[0], GG_ansig.shape[1], GG_ansig.shape[2]))  

    for i in range(GG.shape[1]):
        for j in range(GG.shape[2]):
            EE[:,i,j]=2*np.real(np.fft.ifftshift(np.fft.ifft(np.fft.ifftshift(GG_ansig[:,i,j])))) # get E(t) according analytical signal theory

    return EE
