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

numin=0
numax=300

nu=nu[numin:numax]

traj = np.zeros(shape=(128,128), dtype=float)
charges = np.zeros(shape=(30,len(nu)), dtype=float)

for k in range(numin,numax):

    GG1 = GG[k,:,:]
    # GG1 = GG[15,64-20:64+20,64-20:64+20]
    
    pi=math.pi
    
    Nx = 128 #size of the object in pixels
    Ny = 128 #size of the object in pixels
    
    X = 100 #size of the object in mm
    Y = 100 #size of the object in mm   
    
    # -----------------------------------------------------------------------------
    # fig1 = plt.figure(figsize=(2, 2))
    # ax = fig1.add_axes([0,0,1,1])
    # cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
    #                                                 ['xkcd:white','xkcd:orange','xkcd:black'],
    #                                                 256)
    # ax.set_title('|G/G0|',fontsize=12)
    # im=plt.imshow(abs(GG1),aspect='auto',cmap = cmap, interpolation='none')
    # # ax.set_xlim(nu_lim1,nu_lim2)
    # # ax.set_ylim(x_lim1,x_lim2)
    # # plt.clim(0, 2.5)
    # ax.locator_params(axis='x', nbins=5)
    # cax = fig1.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    # plt.colorbar(im, cax=cax, shrink=2)    
    # # plt.close() 
    # # # -----------------------------------------------------------------------------
    # fig2 = plt.figure(figsize=(2, 2))
    # ax = fig2.add_axes([0,0,1,1])
    # quant_steps = 2056
    # cmap=cm.get_cmap('hsv_r',quant_steps)
    # ax.set_title('|G/G0|',fontsize=12)
    # im=plt.imshow(np.angle(GG1),aspect='auto',cmap = cmap, interpolation='none')
    # # ax.set_xlim(nu_lim1,nu_lim2)
    # # ax.set_ylim(x_lim1,x_lim2)
    # plt.clim(-pi, pi)
    # ax.locator_params(axis='x', nbins=5)
    # cax = fig2.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    # plt.colorbar(im, cax=cax, shrink=2)    
    # # plt.close() 
    
    
    # -----------------------------------------------------------------------------
    # find singular by ampl
    
    A = np.power(abs(GG1), 1/15)
    P = np.angle(GG1)
    
    # fig3 = plt.figure(figsize=(2, 2))
    # ax = fig3.add_axes([0,0,1,1])
    # cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
    #                                                 ['xkcd:white','xkcd:orange','xkcd:black'],
    #                                                 256)
    # ax.set_title('|G/G0|',fontsize=12)
    # im=plt.imshow(A,aspect='auto',cmap = cmap, interpolation='none')
    # # ax.set_xlim(nu_lim1,nu_lim2)
    # # ax.set_ylim(x_lim1,x_lim2)
    # # plt.clim(0.8, 1.1)
    # ax.locator_params(axis='x', nbins=5)
    # cax = fig3.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    # plt.colorbar(im, cax=cax, shrink=2)    
    # # plt.close() 
    sizePs = 6
    NNx = 2*sizePs
    NNy = 2*sizePs
    
    XX = X * (NNx/Nx)
    YY = Y * (NNy/Ny)
                   
    N3d=100
    
    ph=np.linspace(-pi,pi,N3d)
    r=XX/5
    
    
    ch=0
    count=0
    tries=Nx*Ny # deep search
    tries=10000
    while ch != 1 and count<tries:
        result = np.where(A == np.amin(A))
        ind = list(zip(result[0], result[1]))
        
        ind = np.asarray(ind, dtype=np.float)
        indx=int(ind[0,0])
        indy=int(ind[0,1])
        
        crop=40
        if indx<crop or indx>Nx-crop or indy<crop or indy>Ny-crop:
            A[indx,indy]=100
            count=count+1
            continue
        
        # if (A[indx,indy]>A[indx+1,indy]) or (A[indx,indy]>A[indx,indy+1]) or (A[indx,indy]>A[indx+1,indy+1]) or (A[indx,indy]>A[indx-1,indy]) or (A[indx,indy]>A[indx,indy-1]) or (A[indx,indy]>A[indx-1,indy-1]) or (A[indx,indy]>A[indx+1,indy-1]) or (A[indx,indy]>A[indx-1,indy+1]):
            
        #     count=count+1
        #     continue                                
        
        Ps = P[indx-sizePs:indx+sizePs,indy-sizePs:indy+sizePs]
        
        # fig4 = plt.figure(figsize=(2, 2))
        # ax = fig4.add_axes([0,0,1,1])
        # cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
        #                                                 ['xkcd:white','xkcd:orange','xkcd:black'],
        #                                                 256)
        # ax.set_title('|G/G0|',fontsize=12)
        # im=plt.imshow(Ps,aspect='auto',cmap = cmap, interpolation='none')
        # # ax.set_xlim(nu_lim1,nu_lim2)
        # # ax.set_ylim(x_lim1,x_lim2)
        # plt.clim(-3.14, 3.14)
        # ax.locator_params(axis='x', nbins=5)
        # cax = fig4.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        # plt.colorbar(im, cax=cax, shrink=2)    
        # # plt.close() 
        # -----------------------------------------------------------------------------
        # find phase rotation

        PPs=np.zeros(shape=(N3d), dtype=float)
        
        for i,phi in enumerate(ph):    
            PPs[i] = Ps[int(r*np.cos(phi)*NNx/XX)+len(Ps)//2, int(r*np.sin(phi)*NNx/XX)+len(Ps)//2]
        
        PPs=np.unwrap(PPs)
        
        charge = (PPs[-1]-PPs[0])/2/pi
                
        # fig7 = plt.figure(figsize=(2, 2))
        # ax = fig7.add_subplot(111, label="1")
        # plt.plot(PPs,'black')
        # plt.xlabel('azim',fontsize=12)
        # plt.ylabel('phi',fontsize=12)
        # ax.set_ylim(-4*pi,4*pi)
        # ax.locator_params(axis='x', nbins=5)
        # -----------------------------------------------------------------------------
        # draw trajectory
        
        if (charge==1 and PPs[0]<PPs[N3d//4]<PPs[N3d//2]<PPs[N3d//2+N3d//4]<PPs[N3d-1]):  
            traj[indx,indy]=1
            charges[:,k-numin]=1
            ch=1
        elif (charge==-1 and PPs[0]>PPs[N3d//4]>PPs[N3d//2]>PPs[N3d//2+N3d//4]>PPs[N3d-1]):
            traj[indx,indy]=-1
            charges[:,k-numin]=-1
            ch=1
        else:A[indx,indy]=100
            
        count=count+1
        
            
fig5 = plt.figure(figsize=(2, 2))
ax = fig5.add_axes([0,0,1,1])
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                ['xkcd:strong blue','xkcd:white','xkcd:red'],
                                                256)
ax.set_title('|G/G0|',fontsize=12)
im=plt.imshow(traj,aspect='auto',cmap = cmap, interpolation='none')
# ax.set_xlim(nu_lim1,nu_lim2)
# ax.set_ylim(x_lim1,x_lim2)
plt.clim(-1, 1)
ax.locator_params(axis='x', nbins=5)
cax = fig5.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im, cax=cax, shrink=2)    
# plt.close() 


z=0
path="D:\\SCIENCE\\2021\\RNF\\test"
prfx=1

visual.traj(traj,X,Y,Nx,Ny,z,typo,path,prfx)
visual.charges(charges,nu,X,Y,Nx,Ny,z,typo,path,prfx)



# fig6 = plt.figure(figsize=(2, 2))
# ax = fig6.add_axes([0,0,1,1])
# cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
#                                                 ['xkcd:strong blue','xkcd:white','xkcd:red'],
#                                                 256)
# ax.set_title('topological charge',fontsize=8)
# plt.xlabel('ν, THz',fontsize=8)
# im=plt.imshow(charges,extent=[nu[0],nu[nu.shape[0]-1],0,29],aspect=0.005,cmap = cmap, interpolation='none')
# plt.clim(-1, 1)
# ax.locator_params(axis='x', nbins=10)
# ax.locator_params(axis='y', nbins=1)
# plt.rc('xtick',labelsize=8)
# ax.set_yticks([])
# # cax = fig6.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
# plt.colorbar(im, cax=cax, shrink=2)
# fig6.set_dpi(120)    
# # plt.close() 


# fig1 = plt.figure(figsize=(2, 2))
# ax = fig1.add_subplot(111, label="1")
# plt.plot(nu,charges,'co', ms=0.5)
# plt.xlabel('nu, THz',fontsize=12)
# plt.ylabel('charge',fontsize=12)
# ax.locator_params(axis='x', nbins=5)

print("--- %s seconds ---" % (time.time() - start_time))



