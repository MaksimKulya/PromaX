# french+our_Q+french
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

pi=math.pi
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

start_time = time.time()

path="D:\\SCIENCE\\2021\\RNF\\Apr\\09.04\\french+our_Q1+french_15plus\\g\\p\\1z"

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
rho=6.6
Gauss2d=Gauss.Gauss(Nx,Ny,X,Y,rho)

# Wavefront G0 initialization 
GG0=G_beam * Gauss2d

x0=64
y0=64
x_lim1=-50;x_lim2=50
y_lim1=-50;y_lim2=50

nu_lim1=0.01;nu_lim2=2.5
t_lim1=-5;t_lim2=5

numin=0;numax=200

tmin=N//2-70;tmax=N//2+70

z=0
typo='trans';prfx=0;y0;clim1=0;clim2=1;norm=0
visual.A_xnu(GG0,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)

x=np.linspace(-X/2,X/2-X/Nx,Nx)
fwhm_1 = fwhm.half_max_x(x,abs(GG0[np.argmax(abs(G)),:,y0]), level=2)
print("FWHM:{:.3f}".format(fwhm_1))

GG0shift=Gshift.Gshift(GG0,nu,c,z,N0)
EE0=backE.backE(GG0shift,N)

typo='trans';prfx=0;clim1=-1;clim2=1;norm=0
visual.E_xt(EE0,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)
# -----------------------------------------------------------------------------
z_ini=180

nu_cr=c*z_ini/X/(X/Nx)
GG0=Prop.Prop(GG0,nu,X,Y,nmedia,c,z_ini,nu_cr)

typo='trans';prfx=0;clim1=0;clim2=1;norm=0
visual.A_xnu(GG0,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z_ini,path,typo,prfx)

GG0shift=Gshift.Gshift(GG0,nu,c,z_ini,N0)
EE0=backE.backE(GG0shift,N)

typo='trans';prfx=0;clim1=-1;clim2=1;norm=0
visual.E_xt(EE0,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z_ini,path,typo,prfx)
# -----------------------------------------------------------------------------
GG0x=GG0
GG0y=np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
# -----------------------------------------------------------------------------
z=0
# 7 qwp
path_qwp = 'D:\\SCIENCE\\2021\\RNF\\plates\\french.txt'
del_b=0
(GG0x_qwp,GG0y_qwp) = QWP.QWP_r(path_qwp,GG0,X,Y,Nx,Ny,nu,GG0x,GG0y,del_b)
# retard_qwp = QWP.efficency(path_qwp,GG0,X,Y,Nx,Ny,nu)        

typo='trans_y';prfx=1;clim1=0;clim2=1;norm=0
visual.PHnu_xy_2gr(GG0x_qwp,GG0y_qwp,nu,X,Nx,x0,y0,nu_lim1,nu_lim2,nu0,G0,norm,c,z_ini,path,typo,prfx)

z_plus=15

nu_cr=c*z_plus/X/(X/Nx)
GG0x_qwp=Prop.Prop(GG0x_qwp,nu,X,Y,nmedia,c,z_plus,nu_cr)
GG0y_qwp=Prop.Prop(GG0y_qwp,nu,X,Y,nmedia,c,z_plus,nu_cr)
# -----------------------------------------------------------------------------
# Qplate
path_qplate = 'D:\\SCIENCE\\2021\\RNF\\plates\\waveplate_half_04-14_12deg.txt'
del_b=12
(GG0x_qplate, GG0y_qplate) = Qplate.Qplate_360(path_qplate,GG0,X,Y,Nx,Ny,nu,GG0x_qwp,GG0y_qwp,del_b)
# retard_qplate = Qplate.efficency(path_qplate,GG0,X,Y,Nx,Ny,nu) 

typo='trans_y';prfx=2;clim1=0;clim2=1;norm=0
visual.PHnu_xy_2gr(GG0x_qplate,GG0y_qplate,nu,X,Nx,x0,y0,nu_lim1,nu_lim2,nu0,G0,norm,c,z_ini,path,typo,prfx)

z_plus=15

nu_cr=c*z_plus/X/(X/Nx)
GG0x_qplate=Prop.Prop(GG0x_qplate,nu,X,Y,nmedia,c,z_plus,nu_cr)
GG0y_qplate=Prop.Prop(GG0y_qplate,nu,X,Y,nmedia,c,z_plus,nu_cr)
# -----------------------------------------------------------------------------
z_ini = z_ini + 2*z_plus

# 7 qwp
path_qwp = 'D:\\SCIENCE\\2021\\RNF\\plates\\french.txt'
del_b=0
(GG0x_qwp2,GG0y_qwp2) = QWP.QWP_r(path_qwp,GG0,X,Y,Nx,Ny,nu,GG0x_qplate,GG0y_qplate,del_b)
# retard_qwp = QWP.efficency(path_qwp,GG0,X,Y,Nx,Ny,nu)        

typo='trans_y';prfx=3;clim1=0;clim2=1;norm=0
visual.A_xnu(GG0y_qwp2,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z_ini,path,typo,prfx)

typo='trans_y';prfx=3;clim1=0;clim2=1;norm=0
visual.Axy_nu(GG0y_qwp2,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z_ini,typo,path,prfx)
visual.PHxy_nu(GG0y_qwp2,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z_ini,typo,path,prfx)

typo='trans_y';prfx=3;clim1=0;clim2=1;norm=0
visual.PHnu_xy_2gr(GG0x_qwp2,GG0y_qwp2,nu,X,Nx,x0,y0,nu_lim1,nu_lim2,nu0,G0,norm,c,z_ini,path,typo,prfx)

GG0x_qwp2shift=Gshift.Gshift(GG0x_qwp2,nu,c,z_ini,N0)
EE0x_qwp2=backE.backE(GG0x_qwp2shift,N)

GG0y_qwp2shift=Gshift.Gshift(GG0y_qwp2,nu,c,z_ini,N0)
EE0y_qwp2=backE.backE(GG0y_qwp2shift,N)

typo='trans_y';prfx=3;clim1=-1;clim2=1;norm=0
visual.E_xt(EE0y_qwp2,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z_ini,path,typo,prfx)

typo='trans_y';prfx=3;clim1=-1;clim2=1;norm=0
visual.Exy_t(EE0y_qwp2,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,z_ini,typo,path,prfx)
# -----------------------------------------------------------------------------
# diffraction
GG0 = GG0y_qwp2

z_gen=105+50
nu_cr=c*z_gen/X/(X/Nx)
GG0=Prop.Prop(GG0,nu,X,Y,nmedia,c,z_gen,nu_cr)

typo='trans_y';prfx=4;clim1=0;clim2=1;norm=0
visual.A_xnu(GG0,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z_gen,path,typo,prfx)

typo='trans_y';prfx=4;clim1=0;clim2=1;norm=0
visual.Axy_nu(GG0,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z_gen,typo,path,prfx)
visual.PHxy_nu(GG0,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z_gen,typo,path,prfx)

GG0shift=Gshift.Gshift(GG0,nu,c,z_ini+z_gen,N0)
EE0=backE.backE(GG0shift,N)

typo='trans_y';prfx=4;clim1=-1;clim2=1;norm=0
visual.E_xt(EE0,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z_gen,path,typo,prfx)
# -----------------------------------------------------------------------------
F=100
parabola = focus.parabola(GG0,nu,X,Y,c,F)

GG0=GG0*parabola
# -----------------------------------------------------------------------------
# Propagation
zz0_1 = np.linspace(0,25,6)
zz0_2 = np.linspace(5,75,15)

zz1_1 = np.linspace(0,50,6)
zz1_2 = np.linspace(10,150,15)

zz2_1 = np.linspace(0,250,6)
zz2_2 = np.linspace(50,750,15)

zz3_1 = np.linspace(0,500,6)
zz3_2 = np.linspace(100,1500,15)

fwhm_2 = np.zeros(shape=(zz3_1.shape[0]+zz3_2.shape[0]), dtype=float)

zzz1_1 = np.linspace(0,50,2)
zzz1_2 = np.linspace(100,150,2)

for i,z in enumerate(zzz1_1):
    
    nu_cr=c*z/X/(X/Nx)
    GG=Prop.Prop(GG0,nu,X,Y,nmedia,c,z,nu_cr)
    GGz=Prop.Propz(GG0,nu,X,Y,nmedia,c,z,nu_cr)
    # -----------------------------------------------------------------------------
    typo='trans_y';prfx=5;clim1=0;clim2=10;norm=0
    visual.A_xnu(GG,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)
    
    typo='trans_y';prfx=5;clim1=0;clim2=1;norm=0
    visual.Axy_nu(GG,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)
    visual.PHxy_nu(GG,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)
    
    typo='long_y';prfx=5;clim1=0;clim2=1;norm=0
    visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)
    # -----------------------------------------------------------------------------
    # backE
    zE = z_ini + z_gen + z
    GGshift=Gshift.Gshift(GG,nu,c,zE,N0)
    EE=backE.backE(GGshift,N)
    
    GGzshift=Gshift.Gshift(GGz,nu,c,zE,N0)
    EEz=backE.backE(GGzshift,N)
    
    typo='trans_y';prfx=5;clim1=-0.5;clim2=0.5;norm=0
    visual.E_xt(EE,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)
    
    typo='long_y';prfx=5;clim1=-0.05;clim2=0.05;norm=0
    visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)
    
    zlast=z
    # fwhm_2[i] = fwhm.half_max_x(x,abs(GG[np.argmax(abs(G)),:,y0]), level=2)

GG0_2 = GG
GG0z_2 = GGz
# GG0_2[:,Nx//2:Nx,:] = 0
# GG0z_2[:,Nx//2:Nx,:] = 0

z=zzz1_1[-1]
typo='trans_y';prfx=6;clim1=0;clim2=10;norm=0
visual.A_xnu(GG0_2,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)

typo='trans_y';prfx=6;clim1=0;clim2=1;norm=0
visual.Axy_nu(GG0_2,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)
visual.PHxy_nu(GG0_2,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)

typo='long_y';prfx=6;clim1=0;clim2=1;norm=0
visual.A_xnu(GG0z_2,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)

zE = z_ini + z_gen + zlast

GG0_2shift=Gshift.Gshift(GG0_2,nu,c,zE,N0)
EE0_2=backE.backE(GG0_2shift,N)

GG0z_2shift=Gshift.Gshift(GG0z_2,nu,c,zE,N0)
EE0z_2=backE.backE(GG0z_2shift,N)

typo='trans_y';prfx=6;clim1=-0.5;clim2=0.5;norm=0
visual.E_xt(EE0_2,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)

typo='long_y';prfx=6;clim1=-0.05;clim2=0.05;norm=0
visual.E_xt(EE0z_2,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)


for i,z in enumerate(zzz1_2):
    
    nu_cr=c*z/X/(X/Nx)
    GG=Prop.Prop(GG0_2,nu,X,Y,nmedia,c,z,nu_cr)
    GGz=Prop.Propz(GG0_2,nu,X,Y,nmedia,c,z,nu_cr)
    # -----------------------------------------------------------------------------
    typo='trans_y';prfx=6;clim1=0;clim2=10;norm=0
    visual.A_xnu(GG,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)
    
    typo='trans_y';prfx=6;clim1=0;clim2=1;norm=0
    visual.Axy_nu(GG,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)
    visual.PHxy_nu(GG,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx)
        
    typo='long_y';prfx=6;clim1=0;clim2=1;norm=0
    visual.A_xnu(GGz,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx)
    
    # -----------------------------------------------------------------------------
    # backE
    zE = z_ini + z_gen + zlast + z
    GGshift=Gshift.Gshift(GG,nu,c,zE,N0)
    EE=backE.backE(GGshift,N)
    
    GGzshift=Gshift.Gshift(GGz,nu,c,zE,N0)
    EEz=backE.backE(GGzshift,N)
    
    typo='trans_y';prfx=6;clim1=-0.5;clim2=0.5;norm=0
    visual.E_xt(EE,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)
    
    typo='long_y';prfx=6;clim1=-0.05;clim2=0.05;norm=0
    visual.E_xt(EEz,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx)
    
    # fwhm_2[i] = fwhm.half_max_x(x,abs(GG[np.argmax(abs(G)),:,y0]), level=2)

# typo='trans';prfx=3;norm=0
# visual.FWHM_z(fwhm_2,norm,c,zz2,path,typo,prfx)

plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
