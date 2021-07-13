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
def parabola(GG,nu,X,Y,c,F):
    
    Nx=GG.shape[1]
    Ny=GG.shape[2]
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)

    lens = np.zeros(shape=(GG.shape[0], GG.shape[1], GG.shape[2]), dtype=nb.types.complex64)

    for k in range(GG.shape[0]): 
        for i in range(GG.shape[1]):
            for j in range(GG.shape[2]):
                lens[k,i,j]=np.exp((1j)*((x[i]**2)+(y[j]**2))*(2*pi*nu[k]/c)/2/F)


    return lens


# @nb.njit
# def parabola2(GG,nu,rho,X,Y,c,F):
    
#     coeff=rho/X
    
#     Nx=int(GG.shape[1]*coeff)
#     Ny=int(GG.shape[2]*coeff)
#     x=np.linspace(-rho/2,rho/2-rho/Nx,Nx)
#     y=np.linspace(-rho/2,rho/2-rho/Ny,Ny)

#     # lens = np.zeros(shape=(GG.shape[0], GG.shape[1], GG.shape[2]), dtype=nb.types.complex64)
#     lens = np.zeros(shape=(GG.shape[0], x.shape[0], y.shape[0]), dtype=nb.types.complex64)

#     for k in range(GG.shape[0]): 
#         for i in range(x.shape[0]):
#             for j in range(y.shape[0]):
#                 lens[k,i,j]=np.exp((1j)*((x[i]**2)+(y[j]**2))*(2*pi*nu[k]/c)/2/F);

#     # temp = (GG.shape[1]//2) - (lens2.shape[1]//2)
#     # lens = np.pad(lens2, ((0, 0), (temp, temp), (temp, temp)), mode='constant')

#     return lens


def axicone(X,Y,Nx,Ny,nu,c,n,alpha):
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
    B1 = np.zeros(shape=(Nx,Ny), dtype=complex)
    phi_1 = np.zeros(shape=(Nx,Ny), dtype=complex)
    rho_1 = np.zeros(shape=(Nx,Ny), dtype=complex)
    
    krho=2*pi*nu/c*(n-1)*alpha
    
    for i in range(Nx):
        for j in range(Ny):
            B1[i,j]=x[i]-1j*y[j]
            input_num = complex(x[i], y[j])
            rho, phi = cmath.polar(input_num)
            phi1=1*phi
            phi_1[i,j]=phi1
            rho_1[i,j]=rho
    
    # F = rho_1/(n-1)/alpha
    
    
    axicone = np.zeros(shape=(nu.shape[0],Nx,Ny), dtype=complex)
    
    for k in range(krho.shape[0]):
        axicone[k,:,:]=np.exp(1j*krho[k]*rho_1)
    
    return axicone







