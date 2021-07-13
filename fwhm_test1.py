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
# import AS
# import C
import Prop
import backE
import read_mask
import visual
import Gshift
import focus
import LG_butch


pi=math.pi

start_time = time.time()

path="D://SCIENCE//2021//RNF//Feb//beams1//butch//collimated//nomask//1z"

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

# nmedia = 1 # refractive index of propagation media
nmedia = np.ones(300)
N0=1

def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def half_max_x(x, y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[0], half),
            lin_interp(x, y, zero_crossings_i[1], half)]

(E,t,G,nu) = pulse.THz(N,T,tau) # E(t) from initial pulse parameters
# -----------------------------------------------------------------------------
# Low pass and high pass filtering
nucrop=300

nu=nu[(N//2)+1:(N//2)+1+nucrop]
G=G[(N//2)+1:(N//2)+1+nucrop]

nu0=nu[np.argmax(abs(G))]
G0=np.max(abs(G))

# -----------------------------------------------------------------------------
# Beam forming
G_beam=beam.beam(G,Nx,Ny)

# Gauss 2d mask
rho=6.6
Gauss2d=Gauss.Gauss(Nx,Ny,X,Y,rho)

# Wavefront G0 initialization 
GG0=G_beam * Gauss2d

y0=100

z_ini=180
nu_cr=c*z_ini/X/(X/Nx)
GG0=Prop.Prop(GG0,nu,X,Y,nmedia,c,z_ini,nu_cr)


# # find the two crossing points
# hmx = half_max_x(x,abs(GG0[np.argmax(abs(G)),:,y0]))

# # print the answer
# fwhm = hmx[1] - hmx[0]
# print("FWHM:{:.3f}".format(fwhm))

# # a convincing plot
# half = max(abs(GG0[np.argmax(abs(G)),:,y0]))/2.0

# fig3 = plt.figure(figsize=(2, 2))
# plt.plot(x,abs(GG0[np.argmax(abs(G)),:,y0]))
# plt.plot(hmx, [half, half])


# butch (LG)
G00 = np.max(abs(GG0))
GG0 = LG_butch.LG_butch(X,Y,Nx,Ny,GG0,nu,G00,nu00=1.5)

typo='trans';prfx=1;y0;nu_lim1=0.01;nu_lim2=2.5;x_lim1=-50;x_lim2=50;clim1=0;clim2=1;norm=0
visual.A_xnu(GG0,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,0,path,typo,prfx)


# find the two crossing points
hmx = half_max_x(x,abs(GG0[np.argmax(abs(G)),:,y0]))

# print the answer
fwhm = hmx[1] - hmx[0]
print("FWHM:{:.3f}".format(fwhm))

# a convincing plot
half = max(abs(GG0[np.argmax(abs(G)),:,y0]))/2.0


fig3 = plt.figure(figsize=(2, 2))
plt.plot(x,abs(GG0[np.argmax(abs(G)),:,y0]))
plt.plot(hmx, [half, half])





b = abs(GG0[np.argmax(abs(G)),:,y0])
a = np.arange(len(b))


max_y = max(b)  # Find the maximum y value
ii = [i for i in range(len(b)) if b[i] > max_y/2]

fwhm2 = x[max(ii)]-x[min(ii)]
print("FWHM:{:.3f}".format(fwhm2))

orang = np.zeros(shape=(len(x)), dtype=float)
orang[ii[0]:ii[-1]] = max_y/2

fig3 = plt.figure(figsize=(2, 2))
plt.plot(x,abs(GG0[np.argmax(abs(G)),:,y0]))
plt.plot(x, orang)







