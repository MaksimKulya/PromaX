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

path="D:\\SCIENCE\\2021\\RNF\\May\\gen\\KumVsLG\\LG\\ao19\\rho=30\\L1"

# Initial parameters
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

N = 2048 #number of points in pulse temporal form
T = 100  #temporal window in ps
tau = 0.3 #pulse duration in ps

n = 1.46  #refractive index
h = 0  #height of phase object in meters

Nx = 128 #size of the object in pixels
Ny = 128 #size of the object in pixels

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
rho=30
Gauss2d=Gauss.Gauss(Nx,Ny,X,Y,rho)

# Wavefront G0 initialization 
GG0=G_beam * Gauss2d
# -----------------------------------------------------------------------------
x0=64
y0=64
x_lim1=-50;x_lim2=50
y_lim1=-50;y_lim2=50

nu_lim1=0.01;nu_lim2=3
t_lim1=-5;t_lim2=5

numin=0;numax=300

tmin=N//2-70;tmax=N//2+70
# -----------------------------------------------------------------------------
En=np.array([])

z=0; zt=0; prfx=0

nu_cr=c*z/X/(X/Nx)
GG0z=Prop.Propz(GG0,nu,X,Y,nmedia,c,z,nu_cr)

typo='compX';clim1=0;clim2=1;norm=0
visual.A_xnu(GG0,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GG0,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=0;clim2=1;norm=0
visual.A_xnu(GG0z,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GG0z,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# x=np.linspace(-X/2,X/2-X/Nx,Nx)
# fwhm_1 = fwhm.half_max_x(x,abs(GG0[np.argmax(abs(G)),:,y0]), level=2)
# print("FWHM:{:.3f}".format(fwhm_1))

EE0=backE.backE(GG0,N,nu,c,zt,N0)
EE0z=backE.backE(GG0z,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=0
visual.E_xt(EE0,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=-1;clim2=1;norm=0
visual.E_xt(EE0z,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
# -----------------------------------------------------------------------------
z_ini=180; z = z_ini; zt = z_ini; prfx=0

nu_cr=c*z_ini/X/(X/Nx)
GG0=Prop.Prop(GG0,nu,X,Y,nmedia,c,z,nu_cr)
# GG0z=Prop.Propz(GG0,nu,X,Y,nmedia,c,0,0)

# typo='compX';clim1=0;clim2=1;norm=0
# visual.A_xnu(GG0,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
# visual.PH_xnu(GG0,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='compZ';clim1=0;clim2=1;norm=0
# visual.A_xnu(GG0z,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
# visual.PH_xnu(GG0z,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# EE0=backE.backE(GG0,N,nu,c,zt,N0)
# EE0z=backE.backE(GG0z,N,nu,c,zt,N0)

# typo='compX';clim1=-1;clim2=1;norm=0
# visual.E_xt(EE0,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# typo='compZ';clim1=-1;clim2=1;norm=0
# visual.E_xt(EE0z,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
# # -----------------------------------------------------------------------------
GGx = GG0
GGy = GG0

# GGav = np.sqrt(abs(GGx)*abs(GGx) + abs(GGy)*abs(GGy))
# GGsum = GGav.sum()
# En=np.insert(En, len(En),GGsum)

# typo='Gav';clim1=0;clim2=1;norm=0
# visual.A_xnu(GGav,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
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


T1=np.exp(1j*1*phi_1)
B1=abs(B1)
B1=B1/np.max(B1[:])

GGx=GGx*B1*T1
GGy=GGy*B1*T1


prfx=5       

# typo='compX';clim1=0;clim2=1;norm=0
# visual.A_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
# visual.PH_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='compY';clim1=0;clim2=1;norm=0
# visual.A_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
# visual.PH_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='PH_delay';clim1=0;clim2=1;norm=0
# visual.PH_delay(GGx,GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# GGav = np.sqrt(abs(GGx)*abs(GGx) + abs(GGy)*abs(GGy))
# GGsum = GGav.sum()
# En=np.insert(En, len(En),GGsum)

# typo='Gav';clim1=0;clim2=1;norm=0
# visual.A_xnu(GGav,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='Gav';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGav,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# GGz = Prop.Propz_polar(GGx,GGy,nu,X,Y,nmedia,c,0,0)

# typo='compZ';clim1=0;clim2=1;norm=0
# visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
# visual.PH_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='compX';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# typo='compY';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# typo='compZ';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGz,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGz,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# typo='retard2';clim1=0;clim2=1;norm=0
# visual.PHnu_xy_2gr(GGx,GGy,nu,X,Nx,x0,y0,nu_lim1,nu_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# EEx=backE.backE(GGx,N,nu,c,zt,N0)
# EEy=backE.backE(GGy,N,nu,c,zt,N0)
# EEz=backE.backE(GGz,N,nu,c,zt,N0)

# typo='compX';clim1=-1;clim2=1;norm=0
# visual.E_xt(EEx,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# typo='compY';clim1=-1;clim2=1;norm=0
# visual.E_xt(EEy,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# typo='compZ';clim1=-1;clim2=1;norm=0
# visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# I = np.sqrt(EEx*EEx + EEy*EEy)
# typo='Iav';clim1=-1;clim2=1;norm=0
# visual.I_xt(I,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# typo='compX';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EEx,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

# typo='compY';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EEy,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

# typo='compZ';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EEz,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

# typo='Iav';clim1=0;clim2=1;norm=0
# visual.Ixy_t(I,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

# typo='compY'
# numin=0
# numax=300
# deepsearch=Nx*Ny
# cropxy=40
# sizePs=6

# findvortex.fnd(GGy,Nx,Ny,X,Y,nu,numin,numax,deepsearch,cropxy,sizePs,z,path,typo,prfx)
# # # -----------------------------------------------------------------------------
# diffraction
z_gen=105+50-15; z = z_gen; zt = zt + z_gen; prfx=6

nu_cr=c*z/X/(X/Nx)
GGx=Prop.Prop(GGx,nu,X,Y,nmedia,c,z,nu_cr)
GGy=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)
GGz=Prop.Propz_polar(GGx,GGy,nu,X,Y,nmedia,c,0,0)

typo='compX';clim1=0;clim2=1;norm=0
visual.A_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.A_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

GGav = np.sqrt(abs(GGx)*abs(GGx) + abs(GGy)*abs(GGy))
GGsum = GGav.sum()
En=np.insert(En, len(En),GGsum)

typo='Gav';clim1=0;clim2=1;norm=0
visual.A_xnu(GGav,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='Gav';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGav,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

typo='compZ';clim1=0;clim2=1;norm=0
visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='compX';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.Axy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
visual.PHxy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# typo='compZ';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGz,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGz,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

EEx=backE.backE(GGx,N,nu,c,zt,N0)
EEy=backE.backE(GGy,N,nu,c,zt,N0)
EEz=backE.backE(GGz,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=0
visual.E_xt(EEx,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=-1;clim2=1;norm=0
visual.E_xt(EEy,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# typo='compZ';clim1=-1;clim2=1;norm=0
# visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

I = np.sqrt(EEx*EEx + EEy*EEy)
typo='Iav';clim1=-1;clim2=1;norm=0
visual.I_xt(I,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compX';clim1=-0.1;clim2=0.1;norm=0
visual.Exy_t(EEx,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

typo='compY';clim1=-0.1;clim2=0.1;norm=0
visual.Exy_t(EEy,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

# typo='compZ';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EEz,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

typo='Iav';clim1=0;clim2=1;norm=0
visual.Ixy_t(I,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

typo='En'
visual.graph1d(En,path,typo,prfx)

# typo='compY'
# numin=0
# numax=300
# deepsearch=Nx*Ny
# cropxy=40
# sizePs=6

# findvortex.fnd(GGy,Nx,Ny,X,Y,nu,numin,numax,deepsearch,cropxy,sizePs,z,path,typo,prfx)

# typo='compX'
# numin=0
# numax=300
# deepsearch=Nx*Ny
# cropxy=40
# sizePs=6

# findvortex.fnd(GGx,Nx,Ny,X,Y,nu,numin,numax,deepsearch,cropxy,sizePs,z,path,typo,prfx)
# -----------------------------------------------------------------------------
F=500
parabola = focus.parabola(GGy,nu,X,Y,c,F)

GGxs=GGx*parabola
GGys=GGy*parabola

z = F; zt = zt + z; prfx=7

nu_cr=c*z/X/(X/Nx)
GGxs=Prop.Prop(GGxs,nu,X,Y,nmedia,c,z,nu_cr)
GGys=Prop.Prop(GGys,nu,X,Y,nmedia,c,z,nu_cr)
GGzs=Prop.Propz_polar(GGxs,GGys,nu,X,Y,nmedia,c,0,0)

typo='compX';clim1=0;clim2=1;norm=0
visual.A_xnu(GGxs,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGxs,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.A_xnu(GGys,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGys,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=0;clim2=1;norm=0
visual.A_xnu(GGzs,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGzs,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='compX';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGxs,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGxs,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# typo='compY';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGys,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGys,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# typo='compZ';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGzs,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGzs,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

EExs=backE.backE(GGxs,N,nu,c,zt,N0)
EEys=backE.backE(GGys,N,nu,c,zt,N0)
EEzs=backE.backE(GGzs,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=0
visual.E_xt(EExs,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=-1;clim2=1;norm=0
visual.E_xt(EEys,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# typo='compZ';clim1=-1;clim2=1;norm=0
# visual.E_xt(EEzs,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# typo='compX';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EExs,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

# typo='compY';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EEys,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

# typo='compZ';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EEzs,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

zt = zt - z
# -----------------------------------------------------------------------------

a = fwhm.get_fwhm(x,GGy,y0)/2
alpha=a/F/(n-1)

axicone = focus.axicone(X,Y,Nx,Ny,nu,c,n,alpha)

GGxb=GGx*axicone
GGyb=GGy*axicone

z = F; zt = zt + z; prfx=8

nu_cr=c*z/X/(X/Nx)
GGxb=Prop.Prop(GGxb,nu,X,Y,nmedia,c,z,nu_cr)
GGyb=Prop.Prop(GGyb,nu,X,Y,nmedia,c,z,nu_cr)
GGzb=Prop.Propz_polar(GGxb,GGyb,nu,X,Y,nmedia,c,0,0)

typo='compX';clim1=0;clim2=1;norm=0
visual.A_xnu(GGxb,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGxb,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.A_xnu(GGyb,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGyb,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=0;clim2=1;norm=0
visual.A_xnu(GGzb,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGzb,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='compX';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGxb,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGxb,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# typo='compY';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGyb,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGyb,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# typo='compZ';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGzb,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGzb,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

EExb=backE.backE(GGxb,N,nu,c,zt,N0)
EEyb=backE.backE(GGyb,N,nu,c,zt,N0)
EEzb=backE.backE(GGzb,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=0
visual.E_xt(EExb,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=-1;clim2=1;norm=0
visual.E_xt(EEyb,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# typo='compZ';clim1=-1;clim2=1;norm=0
# visual.E_xt(EEzb,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

# typo='compX';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EExb,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

# typo='compY';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EEyb,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

# typo='compZ';clim1=-1;clim2=1;norm=0
# visual.Exy_t(EEzb,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

# sys.modules[__name__].__dict__.clear()

# gc.collect()


