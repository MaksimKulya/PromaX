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

pi=math.pi

Nx = 180 #size of the object in pixels
Ny = 180 #size of the object in pixels

X = 100 #size of the object in mm
Y = 100 #size of the object in mm

x=np.linspace(-X/2,X/2-X/Nx,Nx)
y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)


phi = np.zeros(shape=(Nx,Ny), dtype=float)
rho = np.zeros(shape=(Nx,Ny), dtype=float)

phi15 = np.zeros(shape=(Nx,Ny), dtype=float)
phi9 = np.zeros(shape=(Nx,Ny), dtype=float)
phi4 = np.zeros(shape=(Nx,Ny), dtype=float)
phi8 = np.zeros(shape=(Nx,Ny), dtype=float)
phi3 = np.zeros(shape=(Nx,Ny), dtype=float)
phi7 = np.zeros(shape=(Nx,Ny), dtype=float)
phi2 = np.zeros(shape=(Nx,Ny), dtype=float)
phi6 = np.zeros(shape=(Nx,Ny), dtype=float)


rho1 = np.zeros(shape=(Nx,Ny), dtype=float)
rho2 = np.zeros(shape=(Nx,Ny), dtype=float)

for i in range(Nx):
        for j in range(Ny):
            rho, phi = cmath.polar(complex(x[i], y[j]))
            
            rho1[i,j] = rho
            if 0<=rho<=X/2:
                rho2[i,j]=1
            
            
            if -pi/8<phi<pi/8:
                phi15[i,j]=1
                
            if pi/8<phi<(3*pi/8):
                phi9[i,j]=1
                
            if (3*pi/8)<phi<(5*pi/8):
                phi4[i,j]=1
                
            if (5*pi/8)<phi<(7*pi/8):
                phi8[i,j]=1
                
            if (-7*pi/8)<phi<(7*pi/8):
                phi3[i,j]=1
                                
            if (-7*pi/8)<phi<(-5*pi/8):
                phi7[i,j]=1
                
            if (-5*pi/8)<phi<(-3*pi/8):
                phi2[i,j]=1
                
            if (-3*pi/8)<phi<(-pi/8):
                phi6[i,j]=1


phi3=np.negative(phi3)+1

fig15 = plt.figure(figsize=(2, 2))
ax = fig15.add_axes([0,0,1,1])
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                ['white','xkcd:bright orange','xkcd:chestnut'],
                                                256)
ax.set_title('|G/G0|',fontsize=12)
im=plt.imshow(phi3*rho2,aspect='auto',cmap = cmap, interpolation='none')
# ax.set_xlim(nu_lim1,nu_lim2)
# ax.set_ylim(x_lim1,x_lim2)
# plt.clim(clim1, clim2)
ax.locator_params(axis='x', nbins=5)
cax = fig15.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im, cax=cax, shrink=2)


# fig9 = plt.figure(figsize=(2, 2))
# ax = fig9.add_axes([0,0,1,1])
# cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
#                                                 ['white','xkcd:bright orange','xkcd:chestnut'],
#                                                 256)
# ax.set_title('|G/G0|',fontsize=12)
# im=plt.imshow(phi9*rho2,aspect='auto',cmap = cmap, interpolation='none')
# # ax.set_xlim(nu_lim1,nu_lim2)
# # ax.set_ylim(x_lim1,x_lim2)
# # plt.clim(clim1, clim2)
# ax.locator_params(axis='x', nbins=5)
# cax = fig9.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
# plt.colorbar(im, cax=cax, shrink=2)



# fig8 = plt.figure(figsize=(2, 2))
# ax = fig8.add_axes([0,0,1,1])
# cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
#                                                 ['white','xkcd:bright orange','xkcd:chestnut'],
#                                                 256)
# ax.set_title('|G/G0|',fontsize=12)
# im=plt.imshow(phi3*rho2,aspect='auto',cmap = cmap, interpolation='none')
# # ax.set_xlim(nu_lim1,nu_lim2)
# # ax.set_ylim(x_lim1,x_lim2)
# # plt.clim(clim1, clim2)
# ax.locator_params(axis='x', nbins=5)
# cax = fig8.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
# plt.colorbar(im, cax=cax, shrink=2)










