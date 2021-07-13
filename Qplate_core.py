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
import numba as nb

import Jones

pi=math.pi
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

@nb.njit
def get_jonesmatrix(delta, theta):
    jonesmatrix=np.zeros(shape=(2, 2), dtype=nb.types.complex64)
    jonesmatrix[0,0]=np.cos(delta/2)+1j*np.cos(theta*2)*np.sin(delta/2)
    jonesmatrix[0,1]=1j*np.sin(theta*2)*np.sin(delta/2)
    jonesmatrix[1,0]=1j*np.sin(theta*2)*np.sin(delta/2)
    jonesmatrix[1,1]=np.cos(delta/2)-1j*np.cos(theta*2)*np.sin(delta/2)
    return jonesmatrix


@nb.njit
def loopa(N3d,Nx,Ny,x,y,X,Lf):
    
    sector = np.zeros(shape=(N3d,Nx,Ny), dtype=nb.types.float64)
    sectorx = np.zeros(shape=(N3d,Nx,Ny), dtype=nb.types.float64)
    
    phix = np.zeros(shape=(N3d,Nx,Ny), dtype=nb.types.float64)
    PHI = np.zeros(shape=(N3d,Nx,Ny), dtype=nb.types.float64)
    
    rho1 = np.zeros(shape=(Nx,Ny), dtype=nb.types.float64)
    rho2 = np.zeros(shape=(Nx,Ny), dtype=nb.types.float64)
    
    for k in range(N3d):
        for i in range(Nx):
                for j in range(Ny):
                    rho, phi = cmath.polar(complex(x[i], y[j]))
                    # phi1[i,j] = phi
                    
                    rho1[i,j] = rho
                    if 0<=rho<=X/2:
                        rho2[i,j]=1
                        
                    if -pi <= phi <= -pi+(2*k/(N3d-1))*pi:
                        sector[k,i,j]=1
                        phix[k,i,j]=phi+pi
                          
        sectorx[k,:,:]=sector[k,:,:]-sector[k-1,:,:]  
        PHI[k,:,:]=Lf*(np.max(phix[k,:,:])-np.min(phix[k,:,:]))/2  
        PHI[k,:,:] = PHI[k,:,:] * sectorx[k,:,:]
        
    # take 1st non zero element and put it into new 1d array    
    PHI1 = np.zeros(shape=(N3d), dtype=nb.types.float64)
    for k in range(N3d):    
        for i in range(Nx):
                for j in range(Ny):    
                    if PHI[k,i,j] != 0:
                        PHI1[k] = PHI[k,i,j]
                        break
                    
    return PHI1,sectorx,rho2


@nb.njit
def loopa2(GG0,N3d,nu,h,n_e,n_o,angles,PHI1,del_b,GG0x,GG0y,sectorx,rho2):
    
    GG0x_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=nb.types.complex64)
    GG0y_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=nb.types.complex64) 
    
    GG0x_qplate = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=nb.types.complex64)
    GG0y_qplate = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=nb.types.complex64)
    
    for w in range(N3d):
        
        for k in range(len(nu)):    
            jonesmatrix = np.eye(2, dtype=nb.types.complex64)    
            for q in range(len(h)):     
            
                delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
                jonesmatrice = get_jonesmatrix(delta, angles[q] + PHI1[w] + del_b)
                jonesmatrix= jonesmatrix @ jonesmatrice
                
            GG0x_mod[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*sectorx[w,:,:]*rho2
            GG0y_mod[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*sectorx[w,:,:]*rho2
            
        GG0x_qplate = GG0x_qplate + GG0x_mod 
        GG0y_qplate = GG0y_qplate + GG0y_mod 
            
    return GG0x_qplate,GG0y_qplate

    


