import matplotlib.pyplot as plt
import math
import cmath
import pylab
from matplotlib import mlab
from PIL import Image
import numpy as np
import scipy.constants
import numba as nb

import Cores

pi=math.pi

def AS(GG0,nu,X,Y,nmedia,c,z):

    g = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    GG = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)

    for k in range(GG0.shape[0]):
        g[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(GG0[k,:,:])))  # 2d Fourier to get plane waves spectrum

    (U, H) = Cores.AS_core(g,nu,X,Y,nmedia,c,z)
    
    H[U<0]=0

    gH = g * H

    for k in range(GG0.shape[0]):
        GG[k,:,:] = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(gH[k,:,:])))  # 2d inverse Fourier to get plane waves spectrum

    return GG


def ASz(GG0,nu,X,Y,nmedia,c,z):

    GG0_zp = np.pad(GG0, ((0, 0), (GG0.shape[1]//2, GG0.shape[1]//2), (GG0.shape[1]//2, GG0.shape[1]//2)), mode='constant')
    
    X=X*2
    Y=Y*2    

    g = np.zeros(shape=(GG0_zp.shape[0], GG0_zp.shape[1], GG0_zp.shape[2]), dtype=complex)
    GGz = np.zeros(shape=(GG0_zp.shape[0], GG0_zp.shape[1], GG0_zp.shape[2]), dtype=complex)

    for k in range(GG0.shape[0]):
        g[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(GG0_zp[k,:,:])))  # 2d Fourier to get plane waves spectrum

    
    (U, Hz) = Cores.ASz_core(g,nu,X,Y,nmedia,c,z)
    
    Hz[U<0]=0

    gHz = g * Hz

    for k in range(GG0.shape[0]):
        GGz[k,:,:] = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(gHz[k,:,:])))  # 2d inverse Fourier
        
    GGz = GGz[:,GG0.shape[1]//2:GG0.shape[1]//2+GG0.shape[1],GG0.shape[1]//2:GG0.shape[1]//2+GG0.shape[1]]

    return GGz


def ASz_polar(GG0x,GG0y,nu,X,Y,nmedia,c,z):

    GG0x_zp = np.pad(GG0x, ((0, 0), (GG0x.shape[1]//2, GG0x.shape[1]//2), (GG0x.shape[1]//2, GG0x.shape[1]//2)), mode='constant')
    GG0y_zp = np.pad(GG0y, ((0, 0), (GG0y.shape[1]//2, GG0y.shape[1]//2), (GG0y.shape[1]//2, GG0y.shape[1]//2)), mode='constant')
    
    X=X*2
    Y=Y*2    

    gx = np.zeros(shape=(GG0x_zp.shape[0], GG0x_zp.shape[1], GG0x_zp.shape[2]), dtype=complex)
    gy = np.zeros(shape=(GG0y_zp.shape[0], GG0y_zp.shape[1], GG0y_zp.shape[2]), dtype=complex)
    GGz = np.zeros(shape=(GG0x_zp.shape[0], GG0x_zp.shape[1], GG0x_zp.shape[2]), dtype=complex)

    for k in range(GG0x.shape[0]):
        gx[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(GG0x_zp[k,:,:])))  # 2d Fourier to get plane waves spectrum
        gy[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(GG0y_zp[k,:,:])))  # 2d Fourier to get plane waves spectrum

    (U, Hxz, Hyz) = Cores.ASz_polar_core(gx,gy,nu,X,Y,nmedia,c,z)
    
    Hxz[U<0]=0
    Hyz[U<0]=0

    gHz = gx * Hxz + gy * Hyz

    for k in range(GG0x.shape[0]):
        GGz[k,:,:] = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(gHz[k,:,:])))  # 2d inverse Fourier
        
    GGz = GGz[:,GG0x.shape[1]//2:GG0x.shape[1]//2+GG0x.shape[1],GG0x.shape[1]//2:GG0x.shape[1]//2+GG0x.shape[1]]

    return GGz



def ASeva(GG0,nu,X,Y,nmedia,c,z):

    g = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    GGeva = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)

    for k in range(GG0.shape[0]):
        g[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(GG0[k,:,:])))  # 2d Fourier to get plane waves spectrum

    (U, Heva) = Cores.ASeva_core(g,nu,X,Y,nmedia,c,z)
    
    # Heva[U>0]=0

    gHeva = g * Heva

    for k in range(GG0.shape[0]):
        GGeva[k,:,:] = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(gHeva[k,:,:])))  # 2d inverse Fourier

    return GGeva







