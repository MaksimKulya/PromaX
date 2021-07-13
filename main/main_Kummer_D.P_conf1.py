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

path="D:\\SCIENCE\\2021\\RNF\\May\\diffr_conf1\\B.D.P\\g\\L1\\zz1\\rand_1"

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
rho = 15

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

# -----------------------------------------------------------------------------
z_ini=180; z = z_ini; zt = z_ini; prfx=0

nu_cr=c*z_ini/X/(X/Nx)
GG0=Prop.Prop(GG0,nu,X,Y,nmedia,c,z,nu_cr)

# # -----------------------------------------------------------------------------
GGx = GG0
GGy = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)

GGav = np.sqrt(abs(GGx)*abs(GGx) + abs(GGy)*abs(GGy))
GGsum = GGav.sum()
En=np.insert(En, len(En),GGsum)

# -----------------------------------------------------------------------------
# 7 qwp
path_qwp = 'D:\\SCIENCE\\2021\\RNF\\plates\\new_QWP6.txt'
del_b=0

(GGx,GGy) = QWP.QWP_r(path_qwp,GG0,X,Y,Nx,Ny,nu,GGx,GGy,del_b)   

prfx=1

# -----------------------------------------------------------------------------
z_plus=15; z = z_plus; zt = zt + z_plus; prfx=2

nu_cr=c*z/X/(X/Nx)
GGx=Prop.Prop(GGx,nu,X,Y,nmedia,c,z,nu_cr)
GGy=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)

# -----------------------------------------------------------------------------
# Qplate
path_qplate = 'D:\\SCIENCE\\2021\\RNF\\plates\\waveplate_half_04-14_12deg.txt'
del_b=12

(GGx, GGy) = Qplate.Qplate_360rnb(path_qplate,GG0,X,Y,Nx,Ny,nu,GGx,GGy,del_b,Lf=2)

prfx=3

# -----------------------------------------------------------------------------
z_plus=15; z = z_plus; zt = zt + z_plus; prfx=4

nu_cr=c*z_plus/X/(X/Nx)
GGx=Prop.Prop(GGx,nu,X,Y,nmedia,c,z,nu_cr)
GGy=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)

# # -----------------------------------------------------------------------------
# 7 qwp
path_qwp = 'D:\\SCIENCE\\2021\\RNF\\plates\\new_QWP6.txt'
del_b=0
(GGx,GGy) = QWP.QWP_r(path_qwp,GG0,X,Y,Nx,Ny,nu,GGx,GGy,del_b)

prfx=5       

# -----------------------------------------------------------------------------
# diffraction
z_gen=105+50-15; z = z_gen; zt = zt + z_gen; prfx=6

nu_cr=c*z/X/(X/Nx)
GGx=Prop.Prop(GGx,nu,X,Y,nmedia,c,z,nu_cr)
GGy=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)   
# -----------------------------------------------------------------------------

F=100
parabola = focus.parabola(GGy,nu,X,Y,c,F)

a = fwhm.get_fwhm(x,GGy,y0)/2
alpha=a/F/(n-1)
axicone = focus.axicone(X,Y,Nx,Ny,nu,c,n,alpha)

# GGx=GGx*parabola
# GGy=GGy*parabola

GGx=GGx*axicone
GGy=GGy*axicone
# -----------------------------------------------------------------------------
# Propagation   
zz1_1 = np.linspace(0,F/2,2)
zz1_2 = np.linspace(F/2,3*F/2,3)

for i,z1 in enumerate(zz1_1):
    
    z = z1; zt = zt + z1; prfx=7
    
    nu_cr=c*z/X/(X/Nx)
    GGxf = Prop.Prop(GGx,nu,X,Y,nmedia,c,z,nu_cr)
    GGyf = Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)
    GGzf = Prop.Propz_polar(GGxf,GGyf,nu,X,Y,nmedia,c,0,0)
    # -----------------------------------------------------------------------------
    typo='compX';clim1=0;clim2=10;norm=0
    visual.A_xnu(GGxf,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
    
    typo='compY';clim1=0;clim2=10;norm=0
    visual.A_xnu(GGyf,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
    
    typo='compZ';clim1=0;clim2=1;norm=0
    visual.A_xnu(GGzf,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
    
    
    # typo='compX';clim1=0;clim2=1;norm=0
    # visual.Axy_nu(GGxf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
    # visual.PHxy_nu(GGxf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
    
    # typo='compY';clim1=0;clim2=1;norm=0
    # visual.Axy_nu(GGyf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
    # visual.PHxy_nu(GGyf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
    
    # typo='compZ';clim1=0;clim2=1;norm=0
    # visual.Axy_nu(GGzf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
    # visual.PHxy_nu(GGzf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
 
    EExf=backE.backE(GGxf,N,nu,c,zt,N0)
    EEyf=backE.backE(GGyf,N,nu,c,zt,N0)
    EEzf=backE.backE(GGzf,N,nu,c,zt,N0)
    
    typo='compX';clim1=-0.5;clim2=0.5;norm=0
    visual.E_xt(EExf,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
    
    typo='compY';clim1=-0.5;clim2=0.5;norm=0
    visual.E_xt(EEyf,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
    
    typo='compZ';clim1=-0.05;clim2=0.05;norm=0
    visual.E_xt(EEzf,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
    
    zt = zt - z1

zt = zt + zz1_1[-1]

GGxc = GGxf
GGyc = GGyf
GGzc = GGzf

ampl1 = np.ones(shape=(Nx,Ny), dtype=float)
rand1 = np.random.randint(100, size=(Nx,Ny))/100

n = 1.46  #refractive index
h = 1  #height of phase object in meters
GG_rand1=modulation.modulation(Nx,Ny,ampl1,rand1,nu,n,c,h)

GGxc=GGxc*GG_rand1
GGyc=GGyc*GG_rand1

prfx=8

typo='compX';clim1=0;clim2=10;norm=0
visual.A_xnu(GGxc,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=0;clim2=10;norm=0
visual.A_xnu(GGyc,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=0;clim2=1;norm=0
visual.A_xnu(GGzc,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='compX';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGxc,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGxc,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# typo='compY';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGyc,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGyc,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

# typo='compZ';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGzc,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGzc,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
 
EExc=backE.backE(GGxc,N,nu,c,zt,N0)
EEyc=backE.backE(GGyc,N,nu,c,zt,N0)
EEzc=backE.backE(GGzc,N,nu,c,zt,N0)

typo='compX';clim1=-0.5;clim2=0.5;norm=0
visual.E_xt(EExc,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=-0.5;clim2=0.5;norm=0
visual.E_xt(EEyc,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=-0.05;clim2=0.05;norm=0
visual.E_xt(EEzc,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)


for i,z2 in enumerate(zz1_2):
    
    z = z2; zt = zt + z2; prfx=9
    
    nu_cr=c*z/X/(X/Nx)
    GGxf = Prop.Prop(GGxc,nu,X,Y,nmedia,c,z,nu_cr)
    GGyf = Prop.Prop(GGyc,nu,X,Y,nmedia,c,z,nu_cr)
    GGzf = Prop.Propz_polar(GGxf,GGyf,nu,X,Y,nmedia,c,0,0)
    # -----------------------------------------------------------------------------
    typo='compX';clim1=0;clim2=10;norm=0
    visual.A_xnu(GGxf,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
    
    typo='compY';clim1=0;clim2=10;norm=0
    visual.A_xnu(GGyf,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
    
    typo='compZ';clim1=0;clim2=1;norm=0
    visual.A_xnu(GGzf,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
    
    
    # typo='compX';clim1=0;clim2=1;norm=0
    # visual.Axy_nu(GGxf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
    # visual.PHxy_nu(GGxf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
    
    # typo='compY';clim1=0;clim2=1;norm=0
    # visual.Axy_nu(GGyf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
    # visual.PHxy_nu(GGyf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
    
    # typo='compZ';clim1=0;clim2=1;norm=0
    # visual.Axy_nu(GGzf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
    # visual.PHxy_nu(GGzf,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
 
    EExf=backE.backE(GGxf,N,nu,c,zt,N0)
    EEyf=backE.backE(GGyf,N,nu,c,zt,N0)
    EEzf=backE.backE(GGzf,N,nu,c,zt,N0)
    
    typo='compX';clim1=-0.5;clim2=0.5;norm=0
    visual.E_xt(EExf,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
    
    typo='compY';clim1=-0.5;clim2=0.5;norm=0
    visual.E_xt(EEyf,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
    
    typo='compZ';clim1=-0.05;clim2=0.05;norm=0
    visual.E_xt(EEzf,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
    
    zt = zt - z2

# zt = zt + zz1_2[-1]
    
# typo='compY'
# numin=0
# numax=300
# deepsearch=10000
# cropxy=40
# sizePs=6

# start_time1 = time.time()
# findvortex.fnd(GGyf,Nx,Ny,X,Y,nu,numin,numax,deepsearch,cropxy,sizePs,zt,path,typo,prfx)
# print("--- %s seconds ---" % (time.time() - start_time1))


plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

# sys.modules[__name__].__dict__.clear()

# gc.collect()


