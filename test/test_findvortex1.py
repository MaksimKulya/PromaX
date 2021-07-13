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
from matplotlib import cm

import visual
import findvortex

from numpy import unravel_index

start_time = time.time()

# GG = np.load('D:\\SCIENCE\\2021\\RNF\\Apr\\23.04\\new_QWP6+our_Q+new_QWP6_15plus\\diffr\\rh0=30\\g\\p\\2z\\GG.npy')

GG = np.load('D:\\SCIENCE\\2021\\RNF\\test\\GGy.npy')
nu = np.load('D:\\SCIENCE\\2021\\RNF\\test\\nu.npy')

# nu=np.linspace(1,300,300)
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

typo='compZ';clim1=0;clim2=1;norm=0
# visual.Axy_nu(GG,nu,100,100,128,128,0,300,-50,50,-50,50,1,1,0,c,1,typo,'D:\\SCIENCE\\2021\\RNF\\test2',1)
# visual.PHxy_nu(GG,nu,100,100,128,128,0,300,-50,50,-50,50,1,1,0,c,1,typo,'D:\\SCIENCE\\2021\\RNF\\test2',1)


traj = np.zeros(shape=(128,128), dtype=float)
charges = np.zeros(shape=(30,len(nu)), dtype=float)

typo='compY'
numin=0
numax=300
deepsearch=1000
cropxy=40
sizePs=6

Nx = 128 #size of the object in pixels
Ny = 128 #size of the object in pixels

X = 100 #size of the object in mm
Y = 100 #size of the object in mm  

z=0
path='D:\\SCIENCE\\2021\\RNF\\test'
prfx=1

findvortex.fnd(GG,Nx,Ny,X,Y,nu,numin,numax,deepsearch,cropxy,sizePs,z,path,typo,prfx)

print("--- %s seconds ---" % (time.time() - start_time))












