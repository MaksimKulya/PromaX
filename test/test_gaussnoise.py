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
import findvortex

from scipy.stats import norm
from scipy import signal

np.seterr(divide = 'ignore',invalid='ignore')

pi=math.pi
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

start_time = time.time()

path="D:\\SCIENCE\\2021\\RNF\\May\\diffr_conf1\\vect\\sphere\\g\\L1"

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

noisepow=1/8
(E,t,G,nu) = pulse.THz_noise(N,T,tau,noisepow) # E(t) from initial pulse parameters

# mu=1
# sigma=1

# gauss_fx = norm(mu, sigma)


# noise = np.random.normal(0,1/8,len(E))

# # target_snr_db = 20
# # sig_avg_watts = np.mean(E)
# # sig_avg_db = 10 * np.log10(sig_avg_watts)
# # noise_avg_db = sig_avg_db - target_snr_db
# # noise_avg_watts = 10 ** (noise_avg_db / 10)

# # mean_noise = 0
# # noise_volts = np.random.normal(mean_noise, np.sqrt(noise_avg_watts), len(E))

# E2 = E + noise


fig7 = plt.figure(figsize=(2, 2))
ax = fig7.add_subplot(111, label="1")
plt.plot(E,'black')
plt.xlabel('azim',fontsize=12)
plt.ylabel('phi',fontsize=12)
ax.set_xlim(900,1150)
ax.set_ylim(-1.2,1.2)
ax.locator_params(axis='x', nbins=5)

fig8 = plt.figure(figsize=(2, 2))
ax = fig8.add_subplot(111, label="1")
plt.plot(E2,'black')
plt.xlabel('azim',fontsize=12)
plt.ylabel('phi',fontsize=12)
ax.set_xlim(900,1150)
ax.set_ylim(-1.2,1.2)
ax.locator_params(axis='x', nbins=5)










