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

phi1 = np.zeros(shape=(Nx,Ny), dtype=float)

phi1[10,10]=pi

# for i in range(Nx):
#     for j in range(Ny):
#         rho, phi = cmath.polar(complex(x[i], y[j]))
#         phi1[i,j] = phi

   



fig15 = plt.figure(figsize=(2, 2))
ax = fig15.add_axes([0,0,1,1])
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                ['xkcd:white','xkcd:orange','xkcd:black'],
                                                256)
ax.set_title('|G/G0|',fontsize=12)
im=plt.imshow(phi1,aspect='auto',cmap = cmap, interpolation='none')
# ax.set_xlim(nu_lim1,nu_lim2)
# ax.set_ylim(x_lim1,x_lim2)
plt.clim(0, pi)
ax.locator_params(axis='x', nbins=5)
cax = fig15.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im, cax=cax, shrink=2)    
# plt.close() 













