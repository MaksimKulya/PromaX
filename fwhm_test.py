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

# nmedia = 1 # refractive index of propagation media
nmedia = np.ones(300)
N0=1

(E,t,G,nu) = pulse.THz(N,T,tau) # E(t) from initial pulse parameters

# Low pass and high pass filtering
nucrop=300

nu=nu[(N//2)+1:(N//2)+1+nucrop]
G=G[(N//2)+1:(N//2)+1+nucrop]



# fig2 = plt.figure(figsize=(2, 2))
# ax1=fig2.add_subplot(111, label="1")

# ax1.plot(nu, abs(G),'black')
# ax1.set_xlim(0,3)
# ax1.set_xlabel('ν, THz', fontsize=12, color="black")
# ax1.set_ylabel('|G/G0|', fontsize=12, color="black")
# ax1.tick_params(axis='x', colors="black")
# ax1.tick_params(axis='y', colors="black")
# ax1.locator_params(axis='x', nbins=5)

def peak(x, c):
    return np.exp(-np.power(x - c, 2) / 16.0)

def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def half_max_x(x, y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[0], half),
            lin_interp(x, y, zero_crossings_i[1], half)]



# find the two crossing points
hmx = half_max_x(nu,abs(G))

# print the answer
fwhm = hmx[1] - hmx[0]
print("FWHM:{:.3f}".format(fwhm))

# a convincing plot
half = max(abs(G))/2.0

plt.plot(nu,abs(G))
plt.plot(hmx, [half, half])














