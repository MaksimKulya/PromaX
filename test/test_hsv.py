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

phi1 = np.zeros(shape=(Nx,Ny), dtype=float)

for i in range(Nx):
    for j in range(Ny):
        rho, phi = cmath.polar(complex(x[i], y[j]))
        phi1[i,j] = phi
        
quant_steps = 2056
fig55 = plt.figure(figsize=(3, 3))
ax = fig55.add_subplot(111) 
cmap=cm.get_cmap('hsv_r',quant_steps)
plt.imshow(phi1,aspect='auto',cmap = cmap)
plt.clim(-pi, pi);
# ax.set_figwidth(5)
# ax.set_figheight(4)
#ax.set_aspect(5)
# ax.set_aspect('square')
cax = fig55.add_axes([0.3, 0.1, 0.74, 0.79])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')
# plt.savefig(pt_png+'//'+str(i)+'.png',bbox_inches = 'tight')
# plt.savefig(pt_pdf+'//'+str(i)+'.pdf',bbox_inches = 'tight')
# plt.close()


# fig = plt.figure()

# display_axes = fig.add_axes([0.1,0.1,0.8,0.8], projection='polar')
# display_axes._direction = 2*np.pi ## This is a nasty hack - using the hidden field to 
#                                   ## multiply the values such that 1 become 2*pi
#                                   ## this field is supposed to take values 1 or -1 only!!

# norm = mpl.colors.Normalize(0.0, 2*np.pi)

# # Plot the colorbar onto the polar axis
# # note - use orientation horizontal so that the gradient goes around
# # the wheel rather than centre out
# quant_steps = 2056
# cb = mpl.colorbar.ColorbarBase(display_axes, cmap=cm.get_cmap('hsv',quant_steps),
#                                    norm=norm,
#                                    orientation='horizontal')

# # aesthetics - get rid of border and axis labels                                   
# cb.outline.set_visible(False)                                 
# display_axes.set_axis_off()
# plt.show() # Replace with plt.savefig if you want to save a file

