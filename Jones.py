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
import matplotlib
import matplotlib.pyplot as plt


pi=math.pi
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

# def get_delta(h, n_o, n_e, freq, c):
#     delta = 2*pi*h*1e-3*(n_e-n_o)*freq/c # 1e-3? mm?
#     return delta


def get_jonesmatrix(delta, theta):
    jonesmatrix=np.zeros(shape=(2, 2), dtype=complex)
    jonesmatrix[0,0]=np.cos(delta/2)+1j*np.cos(theta*2)*np.sin(delta/2)
    jonesmatrix[0,1]=1j*np.sin(theta*2)*np.sin(delta/2)
    jonesmatrix[1,0]=1j*np.sin(theta*2)*np.sin(delta/2)
    jonesmatrix[1,1]=np.cos(delta/2)-1j*np.cos(theta*2)*np.sin(delta/2)
    return jonesmatrix


def get_delay(jonesmatrix):
    A=jonesmatrix[0,0]
    B=jonesmatrix[0,1]
    up = A.imag*A.imag+B.imag*B.imag;
    down = A.real*A.real+B.real*B.real
    delay=2*np.arctan(np.sqrt(up/down))
    return delay


def rotation(theta):
    rotation=np.zeros(shape=(2, 2), dtype=complex)
    rotation[0,0]=np.cos(theta)
    rotation[0,1]=-np.sin(theta)
    rotation[1,0]=np.sin(theta)
    rotation[1,1]=np.cos(theta)
    return rotation

# def ExEyJones_Qplate(GG0,GG0x,GG0y,nu,h,n_e,n_o,angles,angle_shift,phi,rho):
    
#     GG0x_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
#     GG0y_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
#     for k in range(len(nu)):    
#         jonesmatrix = np.eye(2, dtype=complex)    
#         for q in range(len(h)):     
        
#             delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
#             jonesmatrice = get_jonesmatrix(delta, angles[q]+angle_shift)
#             jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
#         GG0x_mod[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*phi*rho
#         GG0y_mod[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*phi*rho
        
#         return GG0x_mod,GG0y_mod



# def Mat_mul(Ex, Ey, jonesmatrix):
#     a11=jonesmatrix[0,0]
#     a12=jonesmatrix[0,1]
#     a21=jonesmatrix[1,0]
#     a22=jonesmatrix[1,1]
#     Ex = a11*Ex + a12*Ey
#     Ey = a21*Ex + a22*Ey
#     return Ex,Ey


# def get_retardation(h, angles, n_o, n_e, nu):
#     NN=len(h)
#     retardation=np.zeros(len(nu),1)
#     for ii in range(len(nu)):
#         freq = nu[ii] #*1e12
#         jonesmatrix = np.eye(2, dtype=int)
        
#         for i in range(NN):
#             delta=get_delta(h[i], n_o, n_e, freq)
#             jonesmatrice = get_jonesmatrix(delta, angles[i])
#             jonesmatrix=jonesmatrix*jonesmatrice
#             retardation[ii]=Mat_mul(jonesmatrix)/pi
#     return retardation


# def get_ExEy(Ex,Ey, h, angles, n_o, n_e, nu):
    
#     NN=len(h)
#     retardation=np.zeros(len(nu),1)
    
#     for ii in range(len(nu)):
#         freq = nu[ii] #*1e12
#         jonesmatrix = np.eye(2, dtype=int)
        
#         for i in range(NN):
#             delta=get_delta(h[i], n_o, n_e, freq)
#             jonesmatrice = get_jonesmatrix(delta, angles[i])
#             jonesmatrix=jonesmatrix*jonesmatrice
#             Ex[ii]=Mat_mul(Ex, Ey, jonesmatrix)
#             Ey[ii]=Mat_mul(Ex, Ey, jonesmatrix)
            
#     return Ex, Ey
















