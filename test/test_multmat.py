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

import pulse

pi=math.pi
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



Gx = np.array([1])
Gy = np.array([0])

GGx = np.zeros(shape=(100,100), dtype=float)
GGx[45:55,:]=1
GGy = np.zeros(shape=(100,100), dtype=float)
# GGy[:,45:55]=1


def get_jonesmatrix(delta, theta):
    jonesmatrix=np.zeros(shape=(2, 2), dtype=complex)
    jonesmatrix[0,0]=np.cos(delta/2)+1j*np.cos(theta*2)*np.sin(delta/2)
    jonesmatrix[0,1]=1j*np.sin(theta*2)*np.sin(delta/2)
    jonesmatrix[1,0]=1j*np.sin(theta*2)*np.sin(delta/2)
    jonesmatrix[1,1]=np.cos(delta/2)-1j*np.cos(theta*2)*np.sin(delta/2)
    return jonesmatrix


def rotation(theta):
    rotation=np.zeros(shape=(2, 2), dtype=complex)
    rotation[0,0]=np.cos(theta)
    rotation[0,1]=-np.sin(theta)
    rotation[1,0]=np.sin(theta)
    rotation[1,1]=np.cos(theta)
    return rotation


path_qwp = 'D:\\SCIENCE\\2021\\RNF\\plates\\new_QWP6.txt'

wp = np.loadtxt(path_qwp, comments="#", delimiter=",", unpack=False)

h=wp[:,0]
# h=h[::-1]

angles = wp[:,1]
angles = angles * pi / 180
# angles=angles[::-1]

n_o=2.115
n_e=2.165


nu=1
  
# jonesmatrix = np.eye(2, dtype=complex)   
 
# for q in range(len(h)):     

#     delta = 2*pi*h[q]*(n_e-n_o)*nu/c
#     jonesmatrice = get_jonesmatrix(delta, angles[q])
#     jonesmatrix=jonesmatrix @ jonesmatrice
    
h=3.28
angles=92.677*pi/180

del_b = 30*pi/180

delta = 2*pi*h*(n_e-n_o)*nu/c

jonesmatrice = get_jonesmatrix(delta, angles + del_b)

    
Gx_mod = (jonesmatrice[0,0]*Gx + jonesmatrice[0,1]*Gy)
Gy_mod = (jonesmatrice[1,0]*Gx + jonesmatrice[1,1]*Gy)

G_mod = np.array([Gx_mod,Gy_mod])        


G2 = np.array([Gx,Gy])

Ghat =   rotation(del_b) @ G2 
  
# G=np.sqrt(abs(GGx)**2+abs(GGy)**2)
# G_mod=np.sqrt(abs(GGx_mod)**2+abs(GGy_mod)**2)

# # jonesmatrice2 = get_jonesmatrix(delta, angles)

# jonesmatrice2 = rotation(del_b) @ get_jonesmatrix(delta, angles) @ rotation(-del_b)

    
# Gx_mod2 = (jonesmatrice2[0,0]*Gx + jonesmatrice2[0,1]*Gy)
# Gy_mod2 = (jonesmatrice2[1,0]*Gx + jonesmatrice2[1,1]*Gy)


        






# fig15 = plt.figure(figsize=(2, 2))
# ax = fig15.add_axes([0,0,1,1])
# cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
#                                                 ['xkcd:white','xkcd:orange','xkcd:black'],
#                                                 256)
# ax.set_title('|G/G0|',fontsize=12)
# im=plt.imshow(abs(GGx),aspect='auto',cmap = cmap, interpolation='none')
# # ax.set_xlim(nu_lim1,nu_lim2)
# # ax.set_ylim(x_lim1,x_lim2)
# plt.clim(0, pi)
# ax.locator_params(axis='x', nbins=5)
# cax = fig15.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
# plt.colorbar(im, cax=cax, shrink=2)    
# # plt.close() 


# fig2 = plt.figure(figsize=(2, 2))
# ax = fig2.add_axes([0,0,1,1])
# cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
#                                                 ['xkcd:white','xkcd:orange','xkcd:black'],
#                                                 256)
# ax.set_title('|G/G0|',fontsize=12)
# im=plt.imshow(abs(GGx_mod),aspect='auto',cmap = cmap, interpolation='none')
# # ax.set_xlim(nu_lim1,nu_lim2)
# # ax.set_ylim(x_lim1,x_lim2)
# plt.clim(0, pi)
# ax.locator_params(axis='x', nbins=5)
# cax = fig2.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
# plt.colorbar(im, cax=cax, shrink=2)    
# # plt.close() 








# delta0 = 2*pi*h[0]*(n_e-n_o)*nu/c
# delta1 = 2*pi*h[1]*(n_e-n_o)*nu/c
# delta2 = 2*pi*h[2]*(n_e-n_o)*nu/c
# delta3 = 2*pi*h[3]*(n_e-n_o)*nu/c
# delta4 = 2*pi*h[4]*(n_e-n_o)*nu/c
# delta5 = 2*pi*h[5]*(n_e-n_o)*nu/c
# delta6 = 2*pi*h[6]*(n_e-n_o)*nu/c


# jonesmatrice0 = get_jonesmatrix(delta0, angles[0])
# jonesmatrice1 = get_jonesmatrix(delta1, angles[1])
# jonesmatrice2 = get_jonesmatrix(delta2, angles[2])
# jonesmatrice3 = get_jonesmatrix(delta3, angles[3])
# jonesmatrice4 = get_jonesmatrix(delta4, angles[4])
# jonesmatrice5 = get_jonesmatrix(delta5, angles[5])
# jonesmatrice6 = get_jonesmatrix(delta6, angles[6])


# jonesmatrix2 = jonesmatrice0 @ jonesmatrice1 @ jonesmatrice2 @ jonesmatrice3 @ jonesmatrice4 @ jonesmatrice5 @ jonesmatrice6 

























# a = np.array([[1, 2],[3, 4]], dtype=float)

# b = np.array([[5, 6],[7, 8]], dtype=float)

# c1 = np.matmul(a,b)

# cc1 = a @ b

# c2 = np.matmul(b,a)

# Rinv = np.linalg.inv(R)

# b=R @ b @ np.linalg.inv(R)

# x1 = np.matmul(R,b)
# x1 = np.matmul(x1,Rinv)
