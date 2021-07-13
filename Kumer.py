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

pi=math.pi

def Kumer(X,Y,Nx,Ny,GG0,nu,G0,nu00):
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
    phi_1 = np.zeros(shape=(Nx,Ny), dtype=complex)
    
    for i in range(Nx):
        for j in range(Ny):
            input_num = complex(x[i], y[j])
            rho, phi = cmath.polar(input_num)
            phi1=1*phi
            phi_1[i,j]=phi1
    
    
    T1=np.exp(1j*1*phi_1)
    GG0=GG0*T1
    
    narrowing = np.zeros(shape=(nu.shape[0]), dtype=complex)
    
    # nu00=0.5
    for k in range(nu.shape[0]):
        narrowing[k]=(1/2/np.sqrt(2))*(1-np.exp(1j*pi*nu[k]/nu00))
        GG0[k,:,:]=GG0[k,:,:]*narrowing[k]
    
    GG0 = GG0/np.max(GG0)*G0    
    
    return GG0










