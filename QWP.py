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
from itertools import permutations 
import numba as nb

import Jones

pi=math.pi
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

# qwp plate basic
# qwp plate + reverse
def QWP(path,GG0,X,Y,Nx,Ny,nu,GG0x,GG0y,del_b):
    # 7 q/2
    wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
    h=wp[:,0]
    
    angles = wp[:,1]
    angles = angles * pi / 180
    
    del_b = del_b * pi / 180
    
    n_o=2.115
    n_e=2.165
    
    GG0x_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q])
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
        
        jonesmatrix = Jones.rotation(del_b) @ jonesmatrix @ Jones.rotation(-del_b)
        
        GG0x_mod[k,:,:] = jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:]
        GG0y_mod[k,:,:] = jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:]
    
    return GG0x_mod,GG0y_mod


# qwp plate + reverse
def QWP_r(path,GG0,X,Y,Nx,Ny,nu,GG0x,GG0y,del_b):
    # 7 q/2
    wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
    h=wp[:,0]
    h=h[::-1]
    
    angles = wp[:,1]
    angles = angles * pi / 180
    angles=angles[::-1]
    
    del_b = del_b * pi / 180
    
    n_o=2.115
    n_e=2.165
    
    GG0x_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q])
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
        
        jonesmatrix = Jones.rotation(del_b) @ jonesmatrix @ Jones.rotation(-del_b)
        
        GG0x_mod[k,:,:] = jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:]
        GG0y_mod[k,:,:] = jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:]
    
    return GG0x_mod,GG0y_mod



# # qwp plate basic + reverse + rotation
# def QWP_rr(path,GG0,X,Y,Nx,Ny,nu,GG0x,GG0y,del_b):
#     # 7 q/2
#     wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
#     h=wp[:,0]
#     h=h[::-1]
    
#     angles = wp[:,1]
#     angles = angles * pi / 180
#     angles=angles[::-1]
    
#     n_o=2.115
#     n_e=2.165
    
#     del_b = del_b * pi / 180
    
#     GG0x_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
#     GG0y_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
#     for k in range(len(nu)):    
#         jonesmatrix = np.eye(2, dtype=complex)    
#         for q in range(len(h)):     
        
#             delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
#             jonesmatrice = Jones.get_jonesmatrix(delta, angles[q])
#             jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
        
#         jonesmatrix = Jones.rotation(del_b) @ jonesmatrix @ Jones.rotation(-del_b)
        
#         GG0x_mod[k,:,:] = jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:]
#         GG0y_mod[k,:,:] = jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:]
    
#     return GG0x_mod,GG0y_mod


# qwp under perm loop
def QWP_perm(path,GG0,X,Y,Nx,Ny,nu,GG0x,GG0y,del_b,h,angles):
    # 7 q/2
       
    h=h[::-1]
    angles=angles[::-1]
    
    n_o=2.115
    n_e=2.165
    
    GG0x_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
    
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):    
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q])
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
        
        jonesmatrix = Jones.rotation(del_b*pi/180) @ jonesmatrix @ Jones.rotation(-del_b*pi/180)
        
        GG0x_mod[k,:,:] = jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:]
        GG0y_mod[k,:,:] = jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:]
    
    return GG0x_mod,GG0y_mod


def efficency(path,GG0,X,Y,Nx,Ny,nu):
    # 7 q/2
    wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
    h=wp[:,0]
    
    angles = wp[:,1]
    angles = angles * pi / 180
    
    n_o=2.115
    n_e=2.165
       
    retard = np.zeros(shape=(len(nu),1))
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q])
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
            retard[k] = Jones.get_delay(jonesmatrix)/pi

    
    return retard








def QWP_cw(path,GG0,X,Y,Nx,Ny,nu_cw,GG0x,GG0y,theta):
    # 7 q/2
    wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
    h=wp[:,0]
    
    angles = wp[:,1]
    angles = angles * pi / 180
    
    n_o=2.115
    n_e=2.165
    
    GG0x_mod = np.zeros(shape=(GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod = np.zeros(shape=(GG0.shape[1], GG0.shape[2]), dtype=complex)
    
   
    jonesmatrix = np.eye(2, dtype=complex)    
    for q in range(len(h)):     
    
        delta = 2*pi*h[q]*(n_e-n_o)*nu_cw/c
        jonesmatrice = Jones.get_jonesmatrix(delta, angles[q])
        jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
           
    GG0x_mod = jonesmatrix[0,0]*GG0x + jonesmatrix[0,1]*GG0y
    GG0y_mod = jonesmatrix[1,0]*GG0x + jonesmatrix[1,1]*GG0y
    
    return GG0x_mod,GG0y_mod


def QWP1_1layer_cw(path,GG0,X,Y,Nx,Ny,nu_cw,GG0x,GG0y,theta):
    # 7 q/2
    wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
    h=wp[:,0]
    
    angles = wp[:,1]
    angles = angles * pi / 180
    
    n_o=2.115
    n_e=2.165
    
    GG0x_mod = np.zeros(shape=(GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod = np.zeros(shape=(GG0.shape[1], GG0.shape[2]), dtype=complex)
    
   
    jonesmatrix = np.eye(2, dtype=complex)    
   
    
    delta = 2*pi*h[0]*(n_e-n_o)*nu_cw/c
    jonesmatrice = Jones.get_jonesmatrix(delta, angles[0])
    jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
           
    GG0x_mod = jonesmatrix[0,0]*GG0x + jonesmatrix[0,1]*GG0y
    GG0y_mod = jonesmatrix[1,0]*GG0x + jonesmatrix[1,1]*GG0y
    
    return GG0x_mod,GG0y_mod


def efficency_1layer_cw(path,GG0,X,Y,Nx,Ny,nu_cw):
    # 7 q/2
    wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
    h=wp[:,0]
    
    angles = wp[:,1]
    angles = angles * pi / 180
    
    n_o=2.115
    n_e=2.165
          
    jonesmatrix = np.eye(2, dtype=complex)    
           
    delta = 2*pi*h[0]*(n_e-n_o)*nu_cw/c
    jonesmatrice = Jones.get_jonesmatrix(delta, angles[0])
    jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
    
    delay = Jones.get_delay(jonesmatrix)/pi
    
    # A=jonesmatrix[0,0]
    # B=jonesmatrix[0,1]
    # up = A.imag*A.imag+B.imag*B.imag;
    # down = A.real*A.real+B.real*B.real
    # delay=2*np.arctan(np.sqrt(up/down))/pi
    
    return delay



