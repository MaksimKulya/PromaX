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


# phi1 = np.zeros(shape=(Nx,Ny), dtype=float)
# phi2 = np.zeros(shape=(Nx,Ny), dtype=float)


rho1 = np.zeros(shape=(Nx,Ny), dtype=float)
rho2 = np.zeros(shape=(Nx,Ny), dtype=float)

N3d=100

sector = np.zeros(shape=(N3d,Nx,Ny), dtype=float)
sectorx = np.zeros(shape=(N3d,Nx,Ny), dtype=float)

phix = np.zeros(shape=(N3d,Nx,Ny), dtype=float)
PHI = np.zeros(shape=(N3d,Nx,Ny), dtype=float)

for k in range(N3d):
    for i in range(Nx):
            for j in range(Ny):
                rho, phi = cmath.polar(complex(x[i], y[j]))
                # phi1[i,j] = phi
                
                rho1[i,j] = rho
                if 0<=rho<=X/2:
                    rho2[i,j]=1
                    
                if -pi <= phi <= -pi+(2*k/(N3d-1))*pi:
                    sector[k,i,j]=1
                    phix[k,i,j]=phi+pi
                      
    sectorx[k,:,:]=sector[k,:,:]-sector[k-1,:,:]  
    PHI[k,:,:]=(np.max(phix[k,:,:])-np.min(phix[k,:,:]))/2
    PHI[k,:,:] = PHI[k,:,:] * sectorx[k,:,:]


PHI1 = np.zeros(shape=(N3d), dtype=float)
for k in range(N3d):    
    for i in range(Nx):
            for j in range(Ny):    
                if PHI[k,i,j] != 0:
                    PHI1[k] = PHI[k,i,j]
                    break
    
 
                
    # fig15 = plt.figure(figsize=(2, 2))
    # ax = fig15.add_axes([0,0,1,1])
    # cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
    #                                                 ['white','xkcd:bright orange','xkcd:chestnut'],
    #                                                 256)
    # ax.set_title('|G/G0|',fontsize=12)
    # im=plt.imshow(phi3d[k,:,:],aspect='auto',cmap = cmap, interpolation='none')
    # # ax.set_xlim(nu_lim1,nu_lim2)
    # # ax.set_ylim(x_lim1,x_lim2)
    # plt.clim(0, 1)
    # ax.locator_params(axis='x', nbins=5)
    # cax = fig15.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    # plt.colorbar(im, cax=cax, shrink=2)    
    # # plt.close() 
    
    
    fig15 = plt.figure(figsize=(2, 2))
    ax = fig15.add_axes([0,0,1,1])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                    ['xkcd:white','xkcd:orange','xkcd:black'],
                                                    256)
    ax.set_title('|G/G0|',fontsize=12)
    im=plt.imshow(PHI[k,:,:],aspect='auto',cmap = cmap, interpolation='none')
    # ax.set_xlim(nu_lim1,nu_lim2)
    # ax.set_ylim(x_lim1,x_lim2)
    plt.clim(0, 2*pi)
    ax.locator_params(axis='x', nbins=5)
    cax = fig15.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax, shrink=2)    
    # plt.close() 


# fig15 = plt.figure(figsize=(2, 2))
# ax = fig15.add_axes([0,0,1,1])
# cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
#                                                 ['xkcd:magenta','xkcd:white','xkcd:grassy green'],
#                                                 256)
# ax.set_title('|G/G0|',fontsize=12)
# im=plt.imshow(PHI[10,:,:],aspect='auto',cmap = cmap, interpolation='none')
# # ax.set_xlim(nu_lim1,nu_lim2)
# # ax.set_ylim(x_lim1,x_lim2)
# plt.clim(-pi, pi)
# ax.locator_params(axis='x', nbins=5)
# cax = fig15.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
# plt.colorbar(im, cax=cax, shrink=2)    
# # plt.close() 



del_b=np.linspace(0,2*pi,360)



# for i in range(phi1.shape[0]):
#     for j in range(phi1.shape[1]):
#         if phi1[i,j] == -pi/8:
#             phi2[i,j]=phi1[i,j]


# fig15 = plt.figure(figsize=(2, 2))
# ax = fig15.add_axes([0,0,1,1])
# cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
#                                                 ['white','xkcd:bright orange','xkcd:chestnut'],
#                                                 256)
# ax.set_title('|G/G0|',fontsize=12)
# im=plt.imshow(phi01,aspect='auto',cmap = cmap, interpolation='none')
# # ax.set_xlim(nu_lim1,nu_lim2)
# # ax.set_ylim(x_lim1,x_lim2)
# # plt.clim(clim1, clim2)
# ax.locator_params(axis='x', nbins=5)
# cax = fig15.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
# plt.colorbar(im, cax=cax, shrink=2)


# if -pi/8<phi<pi/8:
#     phi1=1
    


# fig2 = plt.figure(figsize=(2, 2))
# ax = fig2.add_axes([0,0,1,1])
# cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
#                                                 ['white','xkcd:bright orange','xkcd:chestnut'],
#                                                 256)
# ax.set_title('|G/G0|',fontsize=12)
# im=plt.imshow(phi1,aspect='auto',cmap = cmap, interpolation='none')
# # ax.set_xlim(nu_lim1,nu_lim2)
# # ax.set_ylim(x_lim1,x_lim2)
# # plt.clim(clim1, clim2)
# ax.locator_params(axis='x', nbins=5)
# cax = fig2.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
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










