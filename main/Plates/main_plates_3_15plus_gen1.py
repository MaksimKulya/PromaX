# new_QWP6+our_Q+new_QWP6
# diffraction , opaque screen

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

import sys
sys.path[0] = "D://WORK//PromaX"

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

np.seterr(divide = 'ignore',invalid='ignore')

pi=math.pi
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

start_time = time.time()

path="D:\\SCIENCE\\2021\\RNF\\Apr\\generation2\\new_QWP6\\g\\rho=30\\45deg"

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

# nmedia = 1 # refractive index of propagation media
nmedia = np.ones(300)
N0=1

(E,t,G,nu) = pulse.THz(N,T,tau) # E(t) from initial pulse parameters
# -----------------------------------------------------------------------------
# Low pass and high pass filtering
nucrop=300

nu=nu[(N//2)+1:(N//2)+1+nucrop]
G=G[(N//2)+1:(N//2)+1+nucrop]

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

nu_lim1=0.01;nu_lim2=2
t_lim1=-5;t_lim2=5

numin=0;numax=200

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
GG0z=Prop.Propz(GG0,nu,X,Y,nmedia,c,0,0)

typo='compX';clim1=0;clim2=1;norm=0
visual.A_xnu(GG0,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GG0,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=0;clim2=1;norm=0
visual.A_xnu(GG0z,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GG0z,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

EE0=backE.backE(GG0,N,nu,c,zt,N0)
EE0z=backE.backE(GG0z,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=0
visual.E_xt(EE0,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=-1;clim2=1;norm=0
visual.E_xt(EE0z,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
# # -----------------------------------------------------------------------------
GGx = GG0
GGy = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
# -----------------------------------------------------------------------------
# 7 qwp
path_qwp = 'D:\\SCIENCE\\2021\\RNF\\plates\\new_QWP6.txt'
del_b=0
(GGx,GGy) = QWP.QWP_r(path_qwp,GG0,X,Y,Nx,Ny,nu,GGx,GGy,del_b)   

prfx=1

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

typo='PH_delay';clim1=0;clim2=1;norm=0
visual.PH_delay(GGx,GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

nu_cr=c*0/X/(X/Nx)
GGz = Prop.Propz_polar(GGx,GGy,nu,X,Y,nmedia,c,0,0)

typo='compZ';clim1=0;clim2=1;norm=0
visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='retard2';clim1=0;clim2=1;norm=0
visual.PHnu_xy_2gr(GGx,GGy,nu,X,Nx,x0,y0,nu_lim1,nu_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

EEx=backE.backE(GGx,N,nu,c,zt,N0)
EEy=backE.backE(GGy,N,nu,c,zt,N0)
EEz=backE.backE(GGz,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=0
visual.E_xt(EEx,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=-1;clim2=1;norm=0
visual.E_xt(EEy,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=-1;clim2=1;norm=0
visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
# -----------------------------------------------------------------------------
z_plus=15; z = z_plus; zt = zt + z_plus; prfx=2

nu_cr=c*z/X/(X/Nx)
GGx=Prop.Prop(GGx,nu,X,Y,nmedia,c,z,nu_cr)
GGy=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)
GGz = Prop.Propz_polar(GGx,GGy,nu,X,Y,nmedia,c,0,0)

typo='compX';clim1=0;clim2=1;norm=0
visual.A_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.A_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='PH_delay';clim1=0;clim2=1;norm=0
visual.PH_delay(GGx,GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

GGav = np.sqrt(abs(GGx)*abs(GGx) + abs(GGy)*abs(GGy))
GGsum = GGav.sum()
En=np.insert(En, len(En),GGsum)

typo='Gav';clim1=0;clim2=1;norm=0
visual.A_xnu(GGav,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

EEx=backE.backE(GGx,N,nu,c,zt,N0)
EEy=backE.backE(GGy,N,nu,c,zt,N0)
EEz=backE.backE(GGz,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=0
visual.E_xt(EEx,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=-1;clim2=1;norm=0
visual.E_xt(EEy,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=-1;clim2=1;norm=0
visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
# -----------------------------------------------------------------------------
# Qplate
path_qplate = 'D:\\SCIENCE\\2021\\RNF\\plates\\waveplate_half_04-14_12deg.txt'
del_b=12
(GGx, GGy) = Qplate.Qplate_360(path_qplate,GG0,X,Y,Nx,Ny,nu,GGx,GGy,del_b,Lf=2)

prfx=3

typo='compX';clim1=0;clim2=1;norm=0
visual.A_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.A_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compX';clim1=0;clim2=1;norm=0
visual.Axy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
visual.PHxy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.Axy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
visual.PHxy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

GGav = np.sqrt(abs(GGx)*abs(GGx) + abs(GGy)*abs(GGy))
GGsum = GGav.sum()
En=np.insert(En, len(En),GGsum)

typo='Gav';clim1=0;clim2=1;norm=0
visual.A_xnu(GGav,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='PH_delay';clim1=0;clim2=1;norm=0
visual.PH_delay(GGx,GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

GGz = Prop.Propz_polar(GGx,GGy,nu,X,Y,nmedia,c,0,0)

typo='compZ';clim1=0;clim2=1;norm=0
visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='retard2';clim1=0;clim2=1;norm=0
visual.PHnu_xy_2gr(GGx,GGy,nu,X,Nx,x0,y0,nu_lim1,nu_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

EEx=backE.backE(GGx,N,nu,c,zt,N0)
EEy=backE.backE(GGy,N,nu,c,zt,N0)
EEz=backE.backE(GGz,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=0
visual.E_xt(EEx,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=-1;clim2=1;norm=0
visual.E_xt(EEy,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=-1;clim2=1;norm=0
visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
# -----------------------------------------------------------------------------
z_plus=15; z = z_plus; zt = zt + z_plus; prfx=4

nu_cr=c*z_plus/X/(X/Nx)
GGx=Prop.Prop(GGx,nu,X,Y,nmedia,c,z,nu_cr)
GGy=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)
GGz = Prop.Propz_polar(GGx,GGy,nu,X,Y,nmedia,c,0,0)

typo='compX';clim1=0;clim2=1;norm=0
visual.A_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.A_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compX';clim1=0;clim2=1;norm=0
visual.Axy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
visual.PHxy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.Axy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
visual.PHxy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

typo='PH_delay';clim1=0;clim2=1;norm=0
visual.PH_delay(GGx,GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

GGav = np.sqrt(abs(GGx)*abs(GGx) + abs(GGy)*abs(GGy))
GGsum = GGav.sum()
En=np.insert(En, len(En),GGsum)

typo='Gav';clim1=0;clim2=1;norm=0
visual.A_xnu(GGav,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

EEx=backE.backE(GGx,N,nu,c,zt,N0)
EEy=backE.backE(GGy,N,nu,c,zt,N0)
EEz=backE.backE(GGz,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=0
visual.E_xt(EEx,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=-1;clim2=1;norm=0
visual.E_xt(EEy,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=-1;clim2=1;norm=0
visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)
# # -----------------------------------------------------------------------------
# 7 qwp
path_qwp = 'D:\\SCIENCE\\2021\\RNF\\plates\\new_QWP6.txt'
del_b=45
(GGx,GGy) = QWP.QWP_r(path_qwp,GG0,X,Y,Nx,Ny,nu,GGx,GGy,del_b)

prfx=5       

typo='compX';clim1=0;clim2=1;norm=0
visual.A_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.A_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='PH_delay';clim1=0;clim2=1;norm=0
visual.PH_delay(GGx,GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

GGav = np.sqrt(abs(GGx)*abs(GGx) + abs(GGy)*abs(GGy))
GGsum = GGav.sum()
En=np.insert(En, len(En),GGsum)

typo='Gav';clim1=0;clim2=1;norm=0
visual.A_xnu(GGav,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compX';clim1=0;clim2=1;norm=0
visual.Axy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
visual.PHxy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

typo='compY';clim1=0;clim2=1;norm=0
visual.Axy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
visual.PHxy_nu(GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

GGz = Prop.Propz_polar(GGx,GGy,nu,X,Y,nmedia,c,0,0)

typo='compZ';clim1=0;clim2=1;norm=0
visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

# typo='compZ';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GGz,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
# visual.PHxy_nu(GGz,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

typo='retard2';clim1=0;clim2=1;norm=0
visual.PHnu_xy_2gr(GGx,GGy,nu,X,Nx,x0,y0,nu_lim1,nu_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

EEx=backE.backE(GGx,N,nu,c,zt,N0)
EEy=backE.backE(GGy,N,nu,c,zt,N0)
EEz=backE.backE(GGz,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=0
visual.E_xt(EEx,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compY';clim1=-1;clim2=1;norm=0
visual.E_xt(EEy,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=-1;clim2=1;norm=0
visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='compX';clim1=-1;clim2=1;norm=0
visual.Exy_t(EEx,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

typo='compY';clim1=-1;clim2=1;norm=0
visual.Exy_t(EEy,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)

typo='compZ';clim1=-1;clim2=1;norm=0
visual.Exy_t(EEz,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,zt,typo,path,prfx)
# # -----------------------------------------------------------------------------
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

typo='PH_delay';clim1=0;clim2=1;norm=0
visual.PH_delay(GGx,GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

GGav = np.sqrt(abs(GGx)*abs(GGx) + abs(GGy)*abs(GGy))
GGsum = GGav.sum()
En=np.insert(En, len(En),GGsum)

typo='Gav';clim1=0;clim2=1;norm=0
visual.A_xnu(GGav,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compZ';clim1=0;clim2=1;norm=0
visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,zt,path,typo,prfx)
visual.PH_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,zt,path,typo,prfx)

typo='compX';clim1=0;clim2=1;norm=0
visual.Axy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)
visual.PHxy_nu(GGx,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,zt,typo,path,prfx)

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

typo='compZ';clim1=-1;clim2=1;norm=0
visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,zt,path,typo,prfx)

typo='En'
visual.graph1d(En,path,typo,prfx)









# # -----------------------------------------------------------------------------
# F=500
# parabola = focus.parabola(GG0,nu,X,Y,c,F)

# GG0=GG0*parabola
# # -----------------------------------------------------------------------------
# # Propagation
# zzz1_1 = np.linspace(0,50,2)
# zzz1_2 = np.linspace(100,150,2)

# zzz2_1 = np.linspace(0,250,2)
# zzz2_2 = np.linspace(500,750,2)

# zzz3_1 = np.linspace(0,500,2)
# zzz3_2 = np.linspace(1000,1500,2)


# for i,z in enumerate(zzz2_1):
    
#     nu_cr=c*z/X/(X/Nx)
#     GG=Prop.Prop(GG0,nu,X,Y,nmedia,c,z,nu_cr)
#     GGz=Prop.Propz(GG0,nu,X,Y,nmedia,c,z,nu_cr)
#     # -----------------------------------------------------------------------------
#     typo='trans_y';prfx=5;clim1=0;clim2=10;norm=0
#     visual.A_xnu(GG,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)
    
#     typo='trans_y';prfx=5;clim1=0;clim2=1;norm=0
#     visual.Axy_nu(GG,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)
#     visual.PHxy_nu(GG,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)
    
#     typo='long_y';prfx=5;clim1=0;clim2=1;norm=0
#     visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)
#     # -----------------------------------------------------------------------------
#     # backE
#     zE = z_ini + z_gen + z
#     GGshift=Gshift.Gshift(GG,nu,c,zE,N0)
#     EE=backE.backE(GGshift,N)
    
#     GGzshift=Gshift.Gshift(GGz,nu,c,zE,N0)
#     EEz=backE.backE(GGzshift,N)
    
#     typo='trans_y';prfx=5;clim1=-0.5;clim2=0.5;norm=0
#     visual.E_xt(EE,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)
    
#     typo='long_y';prfx=5;clim1=-0.05;clim2=0.05;norm=0
#     visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)
    
#     zlast=z
#     # fwhm_2[i] = fwhm.half_max_x(x,abs(GG[np.argmax(abs(G)),:,y0]), level=2)

# GG0_2 = GG
# GG0z_2 = GGz
# # GG0_2[:,Nx//2:Nx,:] = 0
# # GG0z_2[:,Nx//2:Nx,:] = 0

# z=zzz2_1[-1]

# typo='trans_y';prfx=6;clim1=0;clim2=10;norm=0
# visual.A_xnu(GG0_2,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)

# typo='trans_y';prfx=6;clim1=0;clim2=1;norm=0
# visual.Axy_nu(GG0_2,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)
# visual.PHxy_nu(GG0_2,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)

# typo='long_y';prfx=6;clim1=0;clim2=1;norm=0
# visual.A_xnu(GG0z_2,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)

# zE = z_ini + z_gen + zlast

# GG0_2shift=Gshift.Gshift(GG0_2,nu,c,zE,N0)
# EE0_2=backE.backE(GG0_2shift,N)

# GG0z_2shift=Gshift.Gshift(GG0z_2,nu,c,zE,N0)
# EE0z_2=backE.backE(GG0z_2shift,N)

# typo='trans_y';prfx=6;clim1=-0.5;clim2=0.5;norm=0
# visual.E_xt(EE0_2,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)

# typo='long_y';prfx=6;clim1=-0.05;clim2=0.05;norm=0
# visual.E_xt(EE0z_2,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)


# for i,z in enumerate(zzz2_2):
    
#     nu_cr=c*z/X/(X/Nx)
#     GG=Prop.Prop(GG0_2,nu,X,Y,nmedia,c,z,nu_cr)
#     GGz=Prop.Propz(GG0_2,nu,X,Y,nmedia,c,z,nu_cr)
#     # -----------------------------------------------------------------------------
#     typo='trans_y';prfx=6;clim1=0;clim2=10;norm=0
#     visual.A_xnu(GG,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)
    
#     typo='trans_y';prfx=6;clim1=0;clim2=1;norm=0
#     visual.Axy_nu(GG,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)
#     visual.PHxy_nu(GG,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)
        
#     typo='long_y';prfx=6;clim1=0;clim2=1;norm=0
#     visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)
    
#     # -----------------------------------------------------------------------------
#     # backE
#     zE = z_ini + z_gen + zlast + z
#     GGshift=Gshift.Gshift(GG,nu,c,zE,N0)
#     EE=backE.backE(GGshift,N)
    
#     GGzshift=Gshift.Gshift(GGz,nu,c,zE,N0)
#     EEz=backE.backE(GGzshift,N)
    
#     typo='trans_y';prfx=6;clim1=-0.5;clim2=0.5;norm=0
#     visual.E_xt(EE,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)
    
#     typo='long_y';prfx=6;clim1=-0.05;clim2=0.05;norm=0
#     visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)
    
#     # fwhm_2[i] = fwhm.half_max_x(x,abs(GG[np.argmax(abs(G)),:,y0]), level=2)

# # typo='trans';prfx=3;norm=0
# # visual.FWHM_z(fwhm_2,norm,c,zz2,path,typo,prfx)

plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
