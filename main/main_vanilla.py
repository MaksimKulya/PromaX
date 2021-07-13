import numpy as np
import scipy.constants
import time
import math
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
import AS
import C
import Prop
import backE
import read_mask
import visual
import Gshift

start_time = time.time()

path = "D://SCIENCE//2021//PromaX"  # path for results folder

# Initial parameters
c = scipy.constants.c * 1e-9  # speed of light in mm/ps
pi = math.pi

N = 2048  # number of points in pulse temporal form
T = 100  # temporal window in ps
tau = 0.45  # pulse duration in ps

# n = 1.46  # refractive index
# h = 0  # height of phase object in meters

Nx = 128  # size of the object in pixels
Ny = 128  # size of the object in pixels

X = 102.4  # size of the object in mm
Y = 102.4  # size of the object in mm


(E,t,G,nu) = pulse.THz(N,T,tau)  # E(t) from initial pulse parameters
# -----------------------------------------------------------------------------
# Low pass and high pass filtering
nucrop = 200

nu=nu[(N//2)+1:(N//2)+1+nucrop]
G=G[(N//2)+1:(N//2)+1+nucrop]

nu0=nu[np.argmax(abs(G))]
G0=np.max(abs(G))

t_lim1=-5;t_lim2=5;nu_lim1=0;nu_lim2=2;norm=1
visual.pulse(E,t,G,nu,t_lim1,t_lim2,nu_lim1,nu_lim2,nu0,G0,norm,path)
# -----------------------------------------------------------------------------
# Beam forming
G_beam = beam.beam(G,Nx,Ny)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Gauss 2d mask
rho = 3
Gauss2d = Gauss.Gauss(Nx,Ny,X,Y,rho)
# -----------------------------------------------------------------------------
# Wavefront G0 initialization
GG0 = G_beam * Gauss2d
# G0 = G_beam * G_object * Gauss2d
# -----------------------------------------------------------------------------

GG0x = GG0
GG0y = GG0 * np.exp(1j*pi/2)

# Propagation
z = 33  # propagation distance in mm
zt = z; prfx = 1

nmedia = np.ones(nucrop)
N0 = 1

nu_cr = c*z/X/(X/Nx)
GGx = Prop.Prop(GG0x,nu,X,Y,nmedia,c,z,nu_cr)
GGy = Prop.Prop(GG0y,nu,X,Y,nmedia,c,z,nu_cr)
GGz = Prop.Propz_polar(GGx,GGy,nu,X,Y,nmedia,c,0,0)


# -----------------------------------------------------------------------------
x0=64
y0=64
x_lim1=-50;x_lim2=50
y_lim1=-50;y_lim2=50

nu_lim1=0.01;nu_lim2=3
t_lim1=-5;t_lim2=5

numin=0;numax=300

tmin=N//2-70;tmax=N//2+70

prfx=0

typo='compX';clim1=0;clim2=1;norm=1
visual.A_xnu(GGx,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)

typo='compY';clim1=0;clim2=1;norm=1
visual.A_xnu(GGy,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)

typo='compZ';clim1=0;clim2=1;norm=1
visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)

# -----------------------------------------------------------------------------
# backE
EEx = backE.backE(GGx,N,nu,c,zt,N0)
EEy = backE.backE(GGy,N,nu,c,zt,N0)
EEz = backE.backE(GGz,N,nu,c,zt,N0)

typo='compX';clim1=-1;clim2=1;norm=1
visual.E_xt(EEx,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)

typo='compY';clim1=-1;clim2=1;norm=1
visual.E_xt(EEy,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)

typo='compZ';clim1=-1;clim2=1;norm=1
visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)

print("--- %s seconds ---" % (time.time() - start_time))
