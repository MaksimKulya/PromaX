# new_QWP6+our_Q+new_QWP6
# diffraction , opaque screen

import gc
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

import os

path = os.getcwd()  # Current working directory 
parent = os.path.join(path, os.pardir)  # parent directory 
ROOT_DIR = os.path.abspath(parent)

import sys
sys.path[0] = ROOT_DIR  # path to the root folder

import pulse
import beam
import Gauss
import modulation
import Prop
import backE
import read_mask
import visual
import Gshift
import focus
import LG_butch
import fwhm
import Jones
import Qplate
import QWP
import findvortex

np.seterr(divide = 'ignore',invalid='ignore')

pi=math.pi
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

start_time = time.time()

path="D:\\SCIENCE\\2021\\RNF\\May\\diffr_conf_Yar\\sphere\\ao19\\L3\\noise_0\\c"

# Initial parameters
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

N = 2048 #number of points in pulse temporal form
T = 100  #temporal window in ps
tau = 0.3 #pulse duration in ps

n = 1.46  #refractive index
h = 0  #height of phase object in meters

Nx = 200 #size of the object in pixels
Ny = 200 #size of the object in pixels

X = 100 #size of the object in mm
Y = 100 #size of the object in mm

x=np.linspace(-X/2,X/2-X/Nx,Nx)
y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)

# nmedia = 1 # refractive index of propagation media
nmedia = np.ones(300)
N0=1

(E,t,G,nu) = pulse.THz(N,T,tau) # E(t) from initial pulse parameters
# -----------------------------------------------------------------------------
# Low pass and high pass filtering
nucrop=300

nu=nu[(N//2)+1:(N//2)+1+nucrop]
G=G[(N//2)+1:(N//2)+1+nucrop]

# np.save('D:\\SCIENCE\\2021\\RNF\\test\\nu.npy', nu)

nu0=nu[np.argmax(abs(G))]
G0=np.max(abs(G))

t_lim1=-2;t_lim2=2;nu_lim1=0;nu_lim2=2.5;norm=0
visual.pulse(E,t,G,nu,t_lim1,t_lim2,nu_lim1,nu_lim2,nu0,G0,norm,path)
# -----------------------------------------------------------------------------
# Beam forming
G_beam=beam.beam(G,Nx,Ny)

# Gauss 2d mask
rho = 6.6

Gauss2d=Gauss.Gauss(Nx,Ny,X,Y,rho)

# Wavefront G0 initialization 
GG0=G_beam * Gauss2d

# -----------------------------------------------------------------------------
x0=100
y0=100
x_lim1=-10;x_lim2=10
y_lim1=-10;y_lim2=10

nu_lim1=0.01;nu_lim2=3
t_lim1=-5;t_lim2=5

numin=0;numax=100

tmin=N//2-70;tmax=N//2+70
# -----------------------------------------------------------------------------
En=np.array([])

# -----------------------------------------------------------------------------
z_ini=246-15; z = z_ini; zt = z_ini; prfx=0

nu_cr=c*z_ini/X/(X/Nx)
GG0=Prop.Prop(GG0,nu,X,Y,nmedia,c,z,nu_cr)

# # -----------------------------------------------------------------------------
GGy = GG0
# -----------------------------------------------------------------------------
z_plus=15; z = z_plus; zt = zt + z_plus; prfx=2

nu_cr=c*z/X/(X/Nx)
GGy=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)

# -----------------------------------------------------------------------------
B1 = np.zeros(shape=(Nx,Ny), dtype=complex)
phi_1 = np.zeros(shape=(Nx,Ny), dtype=complex)

for i in range(Nx):
    for j in range(Ny):
        B1[i,j]=x[i]-1j*y[j]
        input_num = complex(x[i], y[j])
        rho, phi = cmath.polar(input_num)
        phi1=1*phi
        phi_1[i,j]=phi1


T1=np.exp(1j*3*phi_1)
B1=abs(B1)
B1=B1/np.max(B1[:])

GGy=GGy*T1

# -----------------------------------------------------------------------------
z_plus=15; z = z_plus; zt = zt + z_plus; prfx=4

nu_cr=c*z_plus/X/(X/Nx)
GGy=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)

# -----------------------------------------------------------------------------
# diffraction
z_gen=60-15; z = z_gen; zt = zt + z_gen; prfx=6

nu_cr=c*z/X/(X/Nx)
GGy=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)   

GGy[:,Nx//2:Nx,:] = 0

z_bladelens=60; z = z_bladelens; zt = zt + z_bladelens; prfx=7

nu_cr=c*z/X/(X/Nx)
GGy=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)  

typo='compY';clim1=0;clim2=10;norm=0
visual.A_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='compY';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# # typo='compZ';clim1=0;clim2=1;norm=0
# # visual.Axy_nu(GGz,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# # visual.PHxy_nu(GGz,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
 
# EEy=backE.backE(GGy,N,nu,c,zt,N0)

# typo='compY';clim1=-0.5;clim2=0.5;norm=0
# visual.E_xt(EEy,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# typo='compY';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EEy,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)
# -----------------------------------------------------------------------------
F=95
parabola = focus.parabola(GGy,nu,X,Y,c,F)

# a = fwhm.get_fwhm(x,GGy,y0)/2
# alpha=a/F/(n-1)
# axicone = focus.axicone(X,Y,Nx,Ny,nu,c,n,alpha)

GGy=GGy*parabola
# -----------------------------------------------------------------------------
# Propagation   
  
z = F; zt = zt + z; prfx=8

nu_cr=c*z/X/(X/Nx)
GGy = Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)
# -----------------------------------------------------------------------------
EEy=backE.backE(GGy,N,nu,c,zt,N0)

noisepow=0

for i in range(GG0.shape[1]):
    for j in range(GG0.shape[2]):
        noise = np.random.normal(0,noisepow,len(E))
        EEy[:,i,j]=EEy[:,i,j]+noise
        
GGy = np.zeros(shape=(N,Nx,Ny), dtype=complex)

for i in range(GG0.shape[1]):
    for j in range(GG0.shape[2]):
        GGy[:,i,j]=np.fft.fftshift(np.fft.fft(np.fft.fftshift(EEy[:,i,j])))

GGy=GGy[(N//2)+1:(N//2)+1+nucrop,:,:]

typo='compY';clim1=0;clim2=10;norm=0
visual.A_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.Axy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
visual.PHxy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# EEy=backE.backE(GGy,N,nu,c,0,N0)

# typo='compY';clim1=-0.5;clim2=0.5;norm=0
# visual.E_xt(EEy,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)


# typo='compY'
# numin=0
# numax=300
# deepsearch=10000
# cropxy=40
# sizePs=6

# start_time1 = time.time()
# findvortex.fnd(GGyf,Nx,Ny,X,Y,nu,numin,numax,deepsearch,cropxy,sizePs,zt,path,typo,prfx)
# print("--- %s seconds ---" % (time.time() - start_time1))


# plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

sys.modules[__name__].__dict__.clear()

# gc.collect()


