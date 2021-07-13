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

def C(GG0,nu,X,Y,nmedia,c,z):

    GG0_zp = np.pad(GG0, ((0, 0), (GG0.shape[1]//2, GG0.shape[1]//2), (GG0.shape[1]//2, GG0.shape[1]//2)), mode='constant')
    
    X=X*2
    Y=Y*2
    
    g = np.zeros(shape=(GG0_zp.shape[0], GG0_zp.shape[1], GG0_zp.shape[2]), dtype=complex)
    GG = np.zeros(shape=(GG0_zp.shape[0], GG0_zp.shape[1], GG0_zp.shape[2]), dtype=complex)
    
    for k in range(GG0_zp.shape[0]):
        g[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(GG0_zp[k,:,:])))  # 2d Fourier to get plane waves spectrum

    h=Cores.C_core(g,nu,X,Y,nmedia,c,z)

    H = np.zeros(shape=(g.shape[0], g.shape[1], g.shape[2]), dtype=complex)
    
    for k in range(g.shape[0]):
        H[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(h[k,:,:])))

    gH = g * H

    for k in range(g.shape[0]):
        GG[k,:,:] = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(gH[k,:,:])))  # 2d inverse Fourier to get plane waves spectrum

    GG = GG[:,GG0.shape[1]//2:GG0.shape[1]//2+GG0.shape[1],GG0.shape[1]//2:GG0.shape[1]//2+GG0.shape[1]]

    return GG



def Cz(GG0,nu,X,Y,nmedia,c,z):

    GG0_zp = np.pad(GG0, ((0, 0), (GG0.shape[1]//2, GG0.shape[1]//2), (GG0.shape[1]//2, GG0.shape[1]//2)), mode='constant')
    
    X=X*2
    Y=Y*2
    
    g = np.zeros(shape=(GG0_zp.shape[0], GG0_zp.shape[1], GG0_zp.shape[2]), dtype=complex)
    GGz = np.zeros(shape=(GG0_zp.shape[0], GG0_zp.shape[1], GG0_zp.shape[2]), dtype=complex)
    
    for k in range(GG0_zp.shape[0]):
        g[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(GG0_zp[k,:,:])))  # 2d Fourier to get plane waves spectrum

    (U,h,kxkz)=Cores.Cz_core(g,nu,X,Y,nmedia,c,z)

    kxkz[U<0]=0

    Hz = np.zeros(shape=(g.shape[0], g.shape[1], g.shape[2]), dtype=complex)
    
    for k in range(g.shape[0]):
        Hz[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(h[k,:,:])))*kxkz[k,:,:]

    gHz = g * Hz

    for k in range(g.shape[0]):
        GGz[k,:,:] = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(gHz[k,:,:])))  # 2d inverse Fourier to get plane waves spectrum

    GGz = GGz[:,GG0.shape[1]//2:GG0.shape[1]//2+GG0.shape[1],GG0.shape[1]//2:GG0.shape[1]//2+GG0.shape[1]]

    return GGz


def Cz_polar(GG0x,GG0y,nu,X,Y,nmedia,c,z):

    GG0x_zp = np.pad(GG0x, ((0, 0), (GG0x.shape[1]//2, GG0x.shape[1]//2), (GG0x.shape[1]//2, GG0x.shape[1]//2)), mode='constant')
    GG0y_zp = np.pad(GG0y, ((0, 0), (GG0y.shape[1]//2, GG0y.shape[1]//2), (GG0y.shape[1]//2, GG0y.shape[1]//2)), mode='constant')
    
    X=X*2
    Y=Y*2
    
    gx = np.zeros(shape=(GG0x_zp.shape[0], GG0x_zp.shape[1], GG0x_zp.shape[2]), dtype=complex)
    gy = np.zeros(shape=(GG0y_zp.shape[0], GG0y_zp.shape[1], GG0y_zp.shape[2]), dtype=complex)
    GGz = np.zeros(shape=(GG0x_zp.shape[0], GG0x_zp.shape[1], GG0x_zp.shape[2]), dtype=complex)
    
    for k in range(GG0x_zp.shape[0]):
        gx[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(GG0x_zp[k,:,:])))  # 2d Fourier to get plane waves spectrum
        gy[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(GG0y_zp[k,:,:])))  # 2d Fourier to get plane waves spectrum

    (U,h,kxkz,kykz)=Cores.Cz_polar_core(gx,gy,nu,X,Y,nmedia,c,z)

    kxkz[U<0]=0
    kykz[U<0]=0

    Hxz = np.zeros(shape=(gx.shape[0], gx.shape[1], gx.shape[2]), dtype=complex)
    Hyz = np.zeros(shape=(gx.shape[0], gx.shape[1], gx.shape[2]), dtype=complex)
    
    for k in range(gx.shape[0]):
        Hxz[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(h[k,:,:])))*kxkz[k,:,:]
        Hyz[k,:,:] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(h[k,:,:])))*kykz[k,:,:]

    gHz = gx * Hxz + gy * Hyz

    for k in range(gx.shape[0]):
        GGz[k,:,:] = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(gHz[k,:,:])))  # 2d inverse Fourier to get plane waves spectrum

    GGz = GGz[:,GG0x.shape[1]//2:GG0x.shape[1]//2+GG0x.shape[1],GG0x.shape[1]//2:GG0x.shape[1]//2+GG0x.shape[1]]

    return GGz




