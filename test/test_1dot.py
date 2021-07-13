import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib as mpl
import math
import cmath

pi=math.pi

Nx = 180 #size of the object in pixels
Ny = 180 #size of the object in pixels

X = 100 #size of the object in mm
Y = 100 #size of the object in mm

x=np.linspace(-X/2,X/2-X/Nx,Nx)
y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)

phi1=np.zeros(shape=(Nx,Ny), dtype=float)
phi1[10:20,10:20]=3
        
quant_steps = 2056

# fig55 = plt.figure(figsize=(3, 3))
# ax = fig55.add_subplot(111) 
# cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
#                                                 ['white','xkcd:bright orange','xkcd:chestnut'],
#                                                 256)
# # cmap=cm.get_cmap('hsv',quant_steps)

# plt.imshow(phi1,extent=[x[0],x[Nx-1],y[Ny-1],y[0]],aspect='auto',cmap = cmap)
# ax.set_xlim(-40,40)
# ax.set_ylim(40,-40)
# # plt.clim(-pi, pi);
# # ax.set_figwidth(5)
# # ax.set_figheight(4)
# # ax.set_aspect(5)
# # ax.set_aspect('square')
# cax = fig55.add_axes([0.3, 0.1, 0.74, 0.79])
# cax.get_xaxis().set_visible(False)
# cax.get_yaxis().set_visible(False)
# cax.patch.set_alpha(0)
# cax.set_frame_on(False)
# plt.colorbar(orientation='vertical')
# fig55.set_dpi(100)
# plt.savefig('D:\\SCIENCE\\2021\\RNF\\Apr\\09.04\\french+our_Q1+french\\1.png',bbox_inches = 'tight')
# # plt.close()


fig55 = plt.figure(figsize=(2, 2))
ax = fig55.add_axes([0,0,1,1])
cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                ['white','xkcd:bright orange','xkcd:chestnut'],
                                                256)

ax.set_title('|G/G0|',fontsize=12)
im=plt.imshow(phi1,extent=[x[0],x[Nx-1],y[Ny-1],y[0]],aspect='auto',cmap = cmap, interpolation='none')
ax.set_xlim(-50,50)
ax.set_ylim(50,-50)
# plt.clim(clim1, clim2)
ax.locator_params(axis='x', nbins=5)
cax = fig55.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im, cax=cax, shrink=2)
fig55.set_dpi(100)
plt.savefig('D:\\SCIENCE\\2021\\RNF\\Apr\\09.04\\french+our_Q1+french\\1.png',bbox_inches = 'tight')
# plt.close()









