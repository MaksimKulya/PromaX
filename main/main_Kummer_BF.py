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

path="D:\\SCIENCE\\2021\\RNF\\May\\focus\\axicone\\\g\\L1\\zz3"

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

zz1 = np.linspace(0,150,16)
zz2 = np.linspace(0,750,16)
zz3 = np.linspace(0,1500,16)

# Gauss 2d mask
rhos = np.array([6.6,10,15,20,25,30])

fwhms_6gr = np.array([])
ksis_6gr = np.array([])

for rh,rho in enumerate(rhos):

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
    F=1000
    
    a = fwhm.get_fwhm(x,GGy,y0)/2
    alpha=a/F/(n-1)
    
    axicone = focus.axicone(X,Y,Nx,Ny,nu,c,n,alpha)
    
    # F=1000
    # parabola = focus.parabola(GGy,nu,X,Y,c,F)
    
    GGy=GGy*axicone
    # -----------------------------------------------------------------------------
    # Propagation   
    fwhms = np.zeros(shape=(zz1.shape[0]), dtype=float)
    ksis = np.zeros(shape=(zz1.shape[0]), dtype=float)
    
    zz=zz3
    
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    
    for i,z in enumerate(zz):
        
        nu_cr=c*z/X/(X/Nx)
        GGyf=Prop.Prop(GGy,nu,X,Y,nmedia,c,z,nu_cr)

        # -----------------------------------------------------------------------------
        path_rho = path + '\\' + str(rho)
        typo='compY';clim1=0;clim2=10;norm=0
        
        visual.A_xnu(GGyf,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path_rho,typo,prfx)
        
        fwhms[i] = fwhm.get_fwhm(x, GGyf, y0)
        
        size_det=18 
        ksis[i] = fwhm.En_ondet(X, Nx, size_det, GGyf, GGy)
        
    fwhms_6gr = np.append(fwhms_6gr,fwhms)
    ksis_6gr = np.append(ksis_6gr,ksis)
    
fwhms_6gr.shape = (fwhms_6gr.size//len(zz), len(zz))
ksis_6gr.shape = (ksis_6gr.size//len(zz), len(zz))

# fwhm_2sum_6gr_rsh = fwhm_2sum_6gr.reshape(6,fwhm_2sum_6gr.shape[0]//6)

typo='FWHM';prfx=6;norm=0
visual.FWHM_z_6gr(fwhms_6gr,norm,c,zz,path,typo,prfx)
np.save(path+'//'+ typo + '//'+'fwhms_6gr.npy', fwhms_6gr)

typo='ksi';prfx=6;norm=0
visual.ksi_z_6gr(ksis_6gr,norm,c,zz,path,typo,prfx)
np.save(path+'//'+ typo + '//'+'ksis_6gr.npy', ksis_6gr)

plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

# sys.modules[__name__].__dict__.clear()

# gc.collect()


