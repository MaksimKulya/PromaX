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
import os
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm

pi=math.pi

def pulse(E,t,G,nu,t_lim1,t_lim2,nu_lim1,nu_lim2,nu0,G0,norm,path):
    
    pt_png=path+'//'+'pulse'+'//'+'png'
    pt_pdf=path+'//'+'pulse'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    T0=1/nu0

    if norm==1:
        nu=nu/nu0
        nu_lim1=nu_lim1/nu0
        nu_lim2=nu_lim2/nu0
        t=t/T0
        t_lim1=t_lim1/T0
        t_lim2=t_lim2/T0
        xlabel_nu='ν/ν0'
        xlabel_t='t/T0'
    else:
        xlabel_nu='ν, THz'
        xlabel_t='t, ps'
        nu=nu
# -----------------------------------------------------------------------------
    fig1 = plt.figure(figsize=(2, 2))
    ax = fig1.add_subplot(111, label="1")
    plt.plot(t, E,'black')
    plt.xlabel(xlabel_t,fontsize=12)
    plt.ylabel('E/E0',fontsize=12)
    ax.set_xlim(t_lim1,t_lim2)
    ax.locator_params(axis='x', nbins=5)
    plt.savefig(pt_png+'//'+'E.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//'+'E.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig1)
    gc.collect()
# -----------------------------------------------------------------------------   
    fig2 = plt.figure(figsize=(2, 2))
    ax1=fig2.add_subplot(111, label="1")
    # ax2=fig2.add_subplot(111, label="2", frame_on=False)

    ax1.plot(nu, abs(G)/G0,'black')
    ax1.set_xlim(nu_lim1,nu_lim2)
    ax1.set_xlabel(xlabel_nu, fontsize=12, color="black")
    ax1.set_ylabel('|G/G0|', fontsize=12, color="black")
    ax1.tick_params(axis='x', colors="black")
    ax1.tick_params(axis='y', colors="black")
    ax1.locator_params(axis='x', nbins=5)

    # ax2.plot(nu, np.angle(G),'b')
    # # ax2.xaxis.tick_top()
    # ax2.yaxis.tick_right()
    # ax2.set_xlim(nu_lim1,nu_lim2)
    # # ax2.set_xlabel(xlabel_nu, color="blue")
    # ax2.set_ylabel('φ, rad', fontsize=12, color="blue")
    # ax2.xaxis.set_label_position('bottom') 
    # ax2.yaxis.set_label_position('right') 
    # ax2.tick_params(axis='x', colors="None")
    # ax2.tick_params(axis='y', colors="blue")
    # ax2.locator_params(axis='x', nbins=5)

    plt.savefig(pt_png+'//'+'A.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//'+'A.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig2)
    gc.collect()

def AMxy(AM,X,Y,Nx,Ny,x_lim1,x_lim2,y_lim1,y_lim2,c,path):
    # Just single 2d image
    
    pt_png=path+'//'+'AMxy'+'//'+'png'
    pt_pdf=path+'//'+'AMxy'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
  
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                    ['white','xkcd:bright orange','xkcd:chestnut'],
                                                    256)
    ax.set_title('AM',fontsize=12)
    plt.xlabel('x,mm',fontsize=12)
    plt.ylabel('y,mm',fontsize=12)
    
    plt.imshow(np.transpose(AM),extent=[x[0],x[Nx-1],y[Ny-1],y[0]],aspect='auto',cmap = cmap)
    ax.set_xlim(x_lim1,x_lim2)
    ax.set_ylim(y_lim1,y_lim2)
    cax = fig.add_axes([0.3, 0.1, 0.74, 0.79])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.colorbar(orientation='vertical')
    plt.savefig(pt_png+'//'+'.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//'+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def A_xnu(GG,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,clim1,clim2,nu0,G0,norm,c,z,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//'+'A_xnu'+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'A_xnu'+'//'+'pdf'
    # pt_txt=path+'//'+typo+'//'+'A_xnu'+'//'+'txt'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)
        
    # if not os.path.isdir(pt_txt):
    #     os.makedirs(pt_txt)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    
    lambda0=c/nu0
    
    if norm==1:
        nu=nu/nu0
        x=x/lambda0
        nu_lim1=nu_lim1/nu0
        nu_lim2=nu_lim2/nu0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        xlabel='ν/ν0'
        ylabel='x/λ0'
    else:
        xlabel='ν, THz'
        ylabel='x,mm'
        nu=nu
        x=x
                
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes([0,0,1,1])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                    ['white','xkcd:bright orange','xkcd:chestnut'],
                                                    256)
    ax.set_title('|G/G0|',fontsize=12)
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel(ylabel,fontsize=12)
    im=plt.imshow(np.transpose(np.abs(GG[:,:,y0])/G0),extent=[nu[0],nu[nu.shape[0]-1],x[Nx-1],x[0]],aspect='auto',
                  cmap = cmap, interpolation='none')
    ax.set_xlim(nu_lim1,nu_lim2)
    ax.set_ylim(x_lim2,x_lim1)
    # plt.clim(clim1, clim2)
    ax.locator_params(axis='x', nbins=5)
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax, shrink=2)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+ str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf+'//'+ str(prfx) + '_' + str(round(z,2))+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()
    
    # np.set_printoptions(threshold=np.inf, linewidth=np.inf)  # turn off summarization, line-wrapping
    # with open(pt_txt + '//' + str(prfx) + '_' + str(round(z,2))+'.txt', 'w') as f:
    #     f.write(np.array2string(np.abs(GG[:,:,y0])/G0, separator=', '))


def Ax_ynu(GG,nu,nu_fix1,nu_fix2,X,Nx,y0,x_lim1,x_lim2,nu0,G0,norm,c,z,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//'+'Ax_ynu'+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'Ax_ynu'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)

    lambda0=c/nu0
    
    if norm==1:
        nu=nu/nu0
        x=x/lambda0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        xlabel='x/λ0'
    else:
        xlabel='x,mm'
        nu=nu
        x=x
        
    fig = plt.figure(figsize=(2, 2))
    ax1=fig.add_subplot(111, label="1")
    ax2=fig.add_subplot(111, label="2", frame_on=False)

    ax1.plot(x, np.abs(GG[nu_fix1,:,y0])/G0,'black')
    ax1.set_xlim(x_lim1,x_lim2)
    ax1.set_xlabel(xlabel, color="black")
    ax1.set_ylabel('|G/G0|', color="black")
    ax1.tick_params(axis='x', colors="black")
    ax1.tick_params(axis='y', colors="black")
    ax1.locator_params(axis='x', nbins=5)

    ax2.plot(x, np.abs(GG[nu_fix2,:,y0])/G0,'red')
    ax2.set_xlim(x_lim1,x_lim2)
    ax2.set_ylabel('|G/G0|', color="None")
    # ax2.set_xlabel(xlabel, color="blue")
    # ax2.tick_params(axis='x', colors="None")
    ax2.tick_params(axis='y', colors="None")
    ax2.locator_params(axis='x', nbins=5)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+ str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf+'//'+ str(prfx) + '_' + str(round(z,2)) + '.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()
    

def Axy_nu(GG,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx):
    
    pt_png=path+'//'+typo+'//'+'Axy_nu'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'Axy_nu'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
    lambda0=c/nu0
    
    if norm==1:
        y=y/lambda0
        x=x/lambda0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        y_lim1=y_lim1/lambda0
        y_lim2=y_lim2/lambda0
        xlabel='x/lambda0'
        ylabel='y/lambda0'
        nunorm=nu/nu0
        nulabel='nu/nu0='
    else:
        xlabel='x,mm'
        ylabel='y,mm'
        nu=nu
        x=x
        nunorm=nu
        nulabel='nu='
    
    for i in range(numin,numax):
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(111)
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                        ['white','xkcd:bright orange','xkcd:chestnut'],
                                                        256)
        ax.set_title(nulabel+str(round((nunorm[i]),2)),fontsize=12)
        plt.xlabel(xlabel,fontsize=12)
        plt.ylabel(ylabel,fontsize=12) 
        plt.imshow(np.power(abs(GG[i,:,:])/G0,1/1),extent=[x[0],x[Nx-1],y[Ny-1],y[0]],aspect='auto',cmap = cmap, interpolation='none')
        ax.set_xlim(x_lim1,x_lim2)
        ax.set_ylim(y_lim2,y_lim1)
        # plt.clim(-pi, pi);
        # ax.set_figwidth(5)
        # ax.set_figheight(4)
        #ax.set_aspect(5)
        # ax.set_aspect('square')
        cax = fig.add_axes([0.3, 0.1, 0.74, 0.79])
        cax.get_xaxis().set_visible(False)
        cax.get_yaxis().set_visible(False)
        cax.patch.set_alpha(0)
        cax.set_frame_on(False)
        plt.colorbar(orientation='vertical')
        fig.set_dpi(80)
        plt.savefig(pt_png+'//'+str(i)+'.png',bbox_inches = 'tight')
        plt.savefig(pt_pdf+'//'+str(i)+'.pdf',bbox_inches = 'tight')
        plt.cla() 
        plt.clf() 
        plt.close('all')   
        plt.close(fig)
        gc.collect()
    
def Anu_xy(GG,nu,X,Nx,x0,y0,nu_lim1,nu_lim2,nu0,G0,norm,c,z,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//'+'Anu_xy'+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'Anu_xy'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    if norm==1:
        nu=nu/nu0
        nu_lim1=nu_lim1/nu0
        nu_lim2=nu_lim2/nu0
        xlabel='ν/ν0'
    else:
        xlabel='ν, THz'
        nu=nu
        
    fig = plt.figure(figsize=(2, 2))
    ax1=fig.add_subplot(111, label="1")
    # ax2=fig2.add_subplot(111, label="2", frame_on=False)

    ax1.plot(nu, np.abs(GG[:,x0,y0])/G0,'black')
    ax1.set_xlim(nu_lim1,nu_lim2)
    ax1.set_xlabel(xlabel, color="black")
    ax1.set_ylabel('|G/G0|', color="black")
    ax1.tick_params(axis='x', colors="black")
    ax1.tick_params(axis='y', colors="black")
    ax1.locator_params(axis='x', nbins=5)

    # ax2.plot(nu, np.angle(GG[:,x0,y0]),'b')
    # # ax2.xaxis.tick_top()
    # ax2.yaxis.tick_right()
    # ax2.set_xlim(nu_lim1,nu_lim2)
    # # ax2.set_xlabel(xlabel, color="blue")
    # ax2.set_ylabel('φ, rad', color="blue")
    # ax2.xaxis.set_label_position('bottom') 
    # ax2.yaxis.set_label_position('right') 
    # ax2.tick_params(axis='x', colors="None")
    # ax2.tick_params(axis='y', colors="blue")
    # ax2.locator_params(axis='x', nbins=5)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+ str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf+'//'+ str(prfx) + '_' + str(round(z,2)) + '.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def Anu_xy_3gr(GG,nu,X,Nx,x0_1,y0_1,x0_2,y0_2,x0_3,y0_3,nu_lim1,nu_lim2,nu0,G0,norm,c,z,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//'+'Anu_xy'+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'Anu_xy'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    if norm==1:
        nu=nu/nu0
        nu_lim1=nu_lim1/nu0
        nu_lim2=nu_lim2/nu0
        xlabel='ν/ν0'
    else:
        xlabel='ν, THz'
        nu=nu
        
    fig = plt.figure(figsize=(2, 2))
    ax1=fig.add_subplot(111, label="1")
    ax2=fig.add_subplot(111, label="2", frame_on=False)
    ax3=fig.add_subplot(111, label="3", frame_on=False)

    ax1.plot(nu, np.abs(GG[:,x0_1,y0_1])/G0,'black')
    ax1.set_xlim(nu_lim1,nu_lim2)
    ax1.set_xlabel(xlabel, color="black")
    ax1.set_ylabel('|G/G0|', color="black")
    ax1.tick_params(axis='x', colors="black")
    ax1.tick_params(axis='y', colors="black")
    ax1.locator_params(axis='x', nbins=5)

    ax2.plot(nu, np.abs(GG[:,x0_2,y0_2]),'b')
    ax2.set_xlim(nu_lim1,nu_lim2)
    ax2.set_ylabel('|G/G0|', color="None")
    # ax2.set_xlabel(xlabel, color="blue")
    # ax2.tick_params(axis='x', colors="None")
    ax2.tick_params(axis='y', colors="None")
    ax2.locator_params(axis='x', nbins=5)
    
    ax3.plot(nu, np.abs(GG[:,x0_3,y0_3]),'r')
    ax3.set_xlim(nu_lim1,nu_lim2)
    ax3.set_ylabel('|G/G0|', color="None")
    # ax2.set_xlabel(xlabel, color="blue")
    # ax2.tick_params(axis='x', colors="None")
    ax3.tick_params(axis='y', colors="None")
    ax3.locator_params(axis='x', nbins=5)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+ str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf+'//'+ str(prfx) + '_' + str(round(z,2)) + '.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def PHnu_xy_2gr(GGx,GGy,nu,X,Nx,x0,y0,nu_lim1,nu_lim2,nu0,G0,norm,c,z,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//'+'PHnu_xy'+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'PHnu_xy'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    if norm==1:
        nu=nu/nu0
        nu_lim1=nu_lim1/nu0
        nu_lim2=nu_lim2/nu0
        xlabel='ν/ν0'
    else:
        xlabel='ν, THz'
        nu=nu
        
    fig = plt.figure(figsize=(2, 2))
    ax1=fig.add_subplot(111, label="1")
    # ax2=fig.add_subplot(111, label="2", frame_on=False)
    ax1.plot(nu,  np.unwrap(np.angle(np.nan_to_num(np.divide(GGy[:,x0,y0],GGx[:,x0,y0]), nan=0, posinf=0, neginf=0))), 'black')
    # ax1.plot(nu,  ((np.angle(GG2[:,x0,y0])-np.angle(GG1[:,x0,y0]))-(pi/1000))  % pi / pi,'black')
    # ax1.plot(nu, np.unwrap( np.angle(GG1[:,x0,y0])-np.angle(GG2[:,x0,y0]) ) % pi /pi,'black')
    ax1.set_xlim(nu_lim1,nu_lim2)
    ax1.set_ylim(-5,20)
    ax1.set_xlabel(xlabel, color="black")
    ax1.set_ylabel('\phi y- \phi x, pi', color="black")
    ax1.tick_params(axis='x', colors="black")
    ax1.tick_params(axis='y', colors="black")
    ax1.locator_params(axis='x', nbins=5)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+ str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf+'//'+ str(prfx) + '_' + str(round(z,2)) + '.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def PH_xnu(GG,nu,X,Nx,y0,nu_lim1,nu_lim2,x_lim1,x_lim2,nu0,G0,norm,c,z,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//'+'PH_xnu'+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'PH_xnu'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    
    lambda0=c/nu0
    
    if norm==1:
        nu=nu/nu0
        x=x/lambda0
        nu_lim1=nu_lim1/nu0
        nu_lim2=nu_lim2/nu0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        xlabel='ν/ν0'
        ylabel='x/λ0'
    else:
        xlabel='ν, THz'
        ylabel='x, mm'
        nu=nu
        x=x
        
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes([0,0,1,1])
    
    quant_steps = 2056
    cmap=cm.get_cmap('hsv_r',quant_steps)
    
    # cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
    #                                                 ['xkcd:magenta','xkcd:white','xkcd:grassy green'],
    #                                                 256)
    ax.set_title('phase, rad',fontsize=12)
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel(ylabel,fontsize=12)
    im=plt.imshow(np.transpose(np.angle(GG[:,:,y0])),extent=[nu[0],nu[nu.shape[0]-1],x[Nx-1],x[0]],aspect='auto',cmap = cmap, interpolation='none')
    ax.set_xlim(nu_lim1,nu_lim2)
    ax.set_ylim(x_lim2,x_lim1)
    plt.clim(-pi, pi)
    ax.locator_params(axis='x', nbins=5)
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax, shrink=2)

    fig.set_dpi(80)
    plt.savefig(pt_png + '//' + str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf + '//' + str(prfx) + '_' + str(round(z,2))+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()
    
    # cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
    #                                                ['xkcd:magenta','xkcd:white','xkcd:grassy green'],
    #                                                256)

    # cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
    #                                                 ['xkcd:white','xkcd:gray','xkcd:black'],
    #                                                 256)

def PHxy_nu(GG,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx):
    
    pt_png=path+'//'+typo+'//'+'PHxy_nu'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'PHxy_nu'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
    lambda0=c/nu0
    
    if norm==1:
        y=y/lambda0
        x=x/lambda0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        y_lim1=y_lim1/lambda0
        y_lim2=y_lim2/lambda0
        xlabel='x/lambda0'
        ylabel='y/lambda0'
        nunorm=nu/nu0
        nulabel='nu/nu0='
    else:
        xlabel='x,mm'
        ylabel='y,mm'
        nu=nu
        x=x
        nunorm=nu
        nulabel='nu='
    
    for i in range(numin,numax):
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(111)
        
        
        # cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
        #                                                 ['mediumorchid', 'blueviolet', 'navy', 'royalblue', 'darkslategrey', 'limegreen', 'darkgreen', 'yellow', 'orange', 'red'],
        #                                                 256)
        # quant_steps = 2056
        
        quant_steps = 2056
        cmap=cm.get_cmap('hsv_r',quant_steps)
        
        ax.set_title(nulabel+str(round((nunorm[i]),2)),fontsize=12)
        plt.xlabel(xlabel,fontsize=12)
        plt.ylabel(ylabel,fontsize=12)   
        plt.imshow(np.angle(GG[i,:,:]),extent=[x[0],x[Nx-1],y[Ny-1],y[0]],aspect='auto',cmap = cmap, interpolation='none')
        ax.set_xlim(x_lim1,x_lim2)
        ax.set_ylim(y_lim2,y_lim1)
        plt.clim(-pi, pi);
        # ax.set_figwidth(5)
        # ax.set_figheight(4)
        #ax.set_aspect(5)
        # ax.set_aspect('square')
        cax = fig.add_axes([0.3, 0.1, 0.74, 0.79])
        cax.get_xaxis().set_visible(False)
        cax.get_yaxis().set_visible(False)
        cax.patch.set_alpha(0)
        cax.set_frame_on(False)
        plt.colorbar(orientation='vertical')
        fig.set_dpi(80)
        plt.savefig(pt_png+'//'+str(i)+'.png',bbox_inches = 'tight')
        plt.savefig(pt_pdf+'//'+str(i)+'.pdf',bbox_inches = 'tight')
        plt.cla() 
        plt.clf() 
        plt.close('all')   
        plt.close(fig)
        gc.collect()


def PH_delay(GGx,GGy,nu,X,Y,Nx,Ny,numin,numax,x_lim1,x_lim2,y_lim1,y_lim2,nu0,G0,norm,c,z,typo,path,prfx):
    
    np.seterr(divide = 'ignore',invalid='ignore')
    
    pt_png=path+'//'+typo+'//'+'PH_delay'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'PH_delay'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
    lambda0=c/nu0
    
    if norm==1:
        y=y/lambda0
        x=x/lambda0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        y_lim1=y_lim1/lambda0
        y_lim2=y_lim2/lambda0
        xlabel='x/lambda0'
        ylabel='y/lambda0'
        nunorm=nu/nu0
        nulabel='nu/nu0='
    else:
        xlabel='x,mm'
        ylabel='y,mm'
        nu=nu
        x=x
        nunorm=nu
        nulabel='nu='
    
    for i in range(numin,numax):
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(111)
        
        
        # cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
        #                                                 ['mediumorchid', 'blueviolet', 'navy', 'royalblue', 'darkslategrey', 'limegreen', 'darkgreen', 'yellow', 'orange', 'red'],
        #                                                 256)
        # quant_steps = 2056
        
        quant_steps = 2056
        cmap=cm.get_cmap('hsv_r',quant_steps)
        
        ax.set_title(nulabel+str(round((nunorm[i]),2)),fontsize=12)
        plt.xlabel(xlabel,fontsize=12)
        plt.ylabel(ylabel,fontsize=12)   
        plt.imshow(np.angle(np.nan_to_num(np.divide(GGy[i,:,:],GGx[i,:,:]), nan=0, posinf=0, neginf=0)),extent=[x[0],x[Nx-1],y[Ny-1],y[0]],aspect='auto',cmap = cmap, interpolation='none')
        ax.set_xlim(x_lim1,x_lim2)
        ax.set_ylim(y_lim2,y_lim1)
        plt.clim(-pi, pi);
        # ax.set_figwidth(5)
        # ax.set_figheight(4)
        #ax.set_aspect(5)
        # ax.set_aspect('square')
        cax = fig.add_axes([0.3, 0.1, 0.74, 0.79])
        cax.get_xaxis().set_visible(False)
        cax.get_yaxis().set_visible(False)
        cax.patch.set_alpha(0)
        cax.set_frame_on(False)
        plt.colorbar(orientation='vertical')
        fig.set_dpi(80)
        plt.savefig(pt_png+'//'+str(i)+'.png',bbox_inches = 'tight')
        plt.savefig(pt_pdf+'//'+str(i)+'.pdf',bbox_inches = 'tight')
        plt.cla() 
        plt.clf() 
        plt.close('all')   
        plt.close(fig)
        gc.collect()


def E_xt(EE,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//'+'E_xt'+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'E_xt'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    
    T0=1/nu0
    lambda0=c/nu0
    
    if norm==1:
        t=t/T0
        x=x/lambda0
        t_lim1=t_lim1/T0
        t_lim2=t_lim2/T0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        xlabel='t/T0'
        ylabel='x/λ0'
    else:
        xlabel='t,ps'
        ylabel='x,mm'
        t=t
        x=x
    
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes([0,0,1,1])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                    ['xkcd:strong blue','white','xkcd:red'],
                                                    256)
    ax.set_title('|E/E0|',fontsize=12)
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel(ylabel,fontsize=12)
    im=plt.imshow(np.transpose(EE[:,:,y0]),extent=[t[0],t[t.shape[0]-1],x[Nx-1],x[0]],aspect='auto',cmap = cmap, interpolation='none')
    ax.set_xlim(t_lim1,t_lim2)
    ax.set_ylim(x_lim2,x_lim1)
    clim=np.max(abs(EE[:,:,y0]))
    plt.clim(-clim, clim);
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax, shrink=2)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+ str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf+'//'+ str(prfx) + '_' + str(round(z,2))+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def I_xt(I,t,X,Nx,y0,t_lim1,t_lim2,x_lim1,x_lim2,clim1,clim2,nu0,norm,c,z,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//'+'I_xt'+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'I_xt'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    
    T0=1/nu0
    lambda0=c/nu0
    
    if norm==1:
        t=t/T0
        x=x/lambda0
        t_lim1=t_lim1/T0
        t_lim2=t_lim2/T0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        xlabel='t/T0'
        ylabel='x/λ0'
    else:
        xlabel='t,ps'
        ylabel='x,mm'
        t=t
        x=x
    
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes([0,0,1,1])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                    ['xkcd:white','xkcd:red'],
                                                    256)
    ax.set_title('|E/E0|',fontsize=12)
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel(ylabel,fontsize=12)
    im=plt.imshow(np.transpose(I[:,:,y0]),extent=[t[0],t[t.shape[0]-1],x[Nx-1],x[0]],aspect='auto',cmap = cmap, interpolation='none')
    ax.set_xlim(t_lim1,t_lim2)
    ax.set_ylim(x_lim2,x_lim1)
    # plt.clim(clim1, clim2);
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax, shrink=2)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+ str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf+'//'+ str(prfx) + '_' + str(round(z,2))+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()    


def Exy_t(EE,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,z,typo,path,prfx):
    
    pt_png=path+'//'+typo+'//'+'Exy_t'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'Exy_t'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
    T0=1/nu0
    lambda0=c/nu0
    
    if norm==1:
        t=t/T0
        x=x/lambda0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        y_lim1=y_lim1/lambda0
        y_lim2=y_lim2/lambda0
        xlabel='x/lambda0'
        ylabel='y/lambda0'
        xlabel='t/T0'
        ylabel='x/λ0'
    else:
        xlabel='x,mm'
        ylabel='y,mm'
        t=t
        x=x
        tnorm=t
        tlabel='t='
    
    for i in range(tmin,tmax):
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(111)
        
        
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                            ['xkcd:strong blue','white','xkcd:red'],
                                                            256)       
        ax.set_title(tlabel+str(round((tnorm[i]),2)),fontsize=12)
        plt.xlabel(xlabel,fontsize=12)
        plt.ylabel(ylabel,fontsize=12)
    
        plt.imshow(EE[i,:,:],extent=[x[0],x[Nx-1],y[Ny-1],y[0]],aspect='auto',cmap = cmap, interpolation='none')
        ax.set_xlim(x_lim1,x_lim2)
        ax.set_ylim(y_lim2,y_lim1)
        plt.clim(clim1, clim2);
        # ax.set_figwidth(5)
        # ax.set_figheight(4)
        #ax.set_aspect(5)
        # ax.set_aspect('square')
        cax = fig.add_axes([0.3, 0.1, 0.74, 0.79])
        cax.get_xaxis().set_visible(False)
        cax.get_yaxis().set_visible(False)
        cax.patch.set_alpha(0)
        cax.set_frame_on(False)
        plt.colorbar(orientation='vertical')
        fig.set_dpi(80)
        plt.savefig(pt_png+'//'+str(i)+'.png',bbox_inches = 'tight')
        plt.savefig(pt_pdf+'//'+str(i)+'.pdf',bbox_inches = 'tight')
        plt.cla() 
        plt.clf() 
        plt.close('all')   
        plt.close(fig)
        gc.collect()


def Ixy_t(I,t,X,Y,Nx,Ny,tmin,tmax,x_lim1,x_lim2,y_lim1,y_lim2,clim1,clim2,nu0,norm,c,z,typo,path,prfx):
    
    pt_png=path+'//'+typo+'//'+'Exy_t'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'Exy_t'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
    T0=1/nu0
    lambda0=c/nu0
    
    if norm==1:
        t=t/T0
        x=x/lambda0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        y_lim1=y_lim1/lambda0
        y_lim2=y_lim2/lambda0
        xlabel='x/lambda0'
        ylabel='y/lambda0'
        xlabel='t/T0'
        ylabel='x/λ0'
    else:
        xlabel='x,mm'
        ylabel='y,mm'
        t=t
        x=x
        tnorm=t
        tlabel='t='
    
    for i in range(tmin,tmax):
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(111)
        
        
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                            ['xkcd:white','xkcd:red'],
                                                            256)       
        ax.set_title(tlabel+str(round((tnorm[i]),2)),fontsize=12)
        plt.xlabel(xlabel,fontsize=12)
        plt.ylabel(ylabel,fontsize=12)
    
        plt.imshow(I[i,:,:],extent=[x[0],x[Nx-1],y[Ny-1],y[0]],aspect='auto',cmap = cmap, interpolation='none')
        ax.set_xlim(x_lim1,x_lim2)
        ax.set_ylim(y_lim2,y_lim1)
        plt.clim(clim1, clim2);
        # ax.set_figwidth(5)
        # ax.set_figheight(4)
        #ax.set_aspect(5)
        # ax.set_aspect('square')
        cax = fig.add_axes([0.3, 0.1, 0.74, 0.79])
        cax.get_xaxis().set_visible(False)
        cax.get_yaxis().set_visible(False)
        cax.patch.set_alpha(0)
        cax.set_frame_on(False)
        plt.colorbar(orientation='vertical')
        fig.set_dpi(80)
        plt.savefig(pt_png+'//'+str(i)+'.png',bbox_inches = 'tight')
        plt.savefig(pt_pdf+'//'+str(i)+'.pdf',bbox_inches = 'tight')
        plt.cla() 
        plt.clf() 
        plt.close('all')   
        plt.close(fig)
        gc.collect()        

def Et_xy(EE,t,X,Nx,x0,y0,t_lim1,t_lim2,clim1,clim2,nu0,norm,c,z,path,typo):
    
    pt_png=path+'//'+typo+'//'+'Et_xy'+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'Et_xy'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    T0=1/nu0
    
    if norm==1:
        t=t/T0
        t_lim1=t_lim1/T0
        t_lim2=t_lim2/T0
        xlabel='ν/ν0'
    else:
        xlabel='ν, THz'
        t=t
        
    fig = plt.figure(figsize=(2, 2))
    ax1=fig.add_subplot(111, label="1")
    # ax2=fig2.add_subplot(111, label="2", frame_on=False)
    ax1.plot(t, EE[:,x0,y0],'black')
    ax1.set_xlim(t_lim1,t_lim2)
    # ax1.set_ylim(clim1, clim2);
    ax1.set_xlabel(xlabel, color="black")
    ax1.set_ylabel('|G/G0|', color="black")
    ax1.tick_params(axis='x', colors="black")
    ax1.tick_params(axis='y', colors="black")
    ax1.locator_params(axis='x', nbins=5)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf+'//'+str(round(z,2))+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def V_rnu(V,nu,rr,nu_lim1,nu_lim2,r_lim1,r_lim2,nu0,G0,norm,c,z,path):
    
    pt_png=path+'//'+'V_rnu'+'//'+'png'
    pt_pdf=path+'//'+'V_rnu'+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)
   
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                    ['white','xkcd:bright orange','xkcd:chestnut'],
                                                    256)
    ax.set_title('|G/G0|',fontsize=12)
    plt.xlabel('nu,THz',fontsize=12)
    plt.ylabel('R,mm',fontsize=12)
    plt.imshow(np.transpose(np.abs(V)/G0),extent=[nu[0],nu[nu.shape[0]-1],rr[rr.shape[0]-1],rr[0]],aspect='auto',cmap = cmap, interpolation='none')
    ax.set_xlim(nu_lim1,nu_lim2)
    ax.set_ylim(r_lim2,r_lim1)
    #plt.clim(-1, 1);
    cax = fig.add_axes([0.3, 0.1, 0.74, 0.79])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.colorbar(orientation='vertical')
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+str(z)+'.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//'+str(z)+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def eta(eta,nu,nu_lim1,nu_lim2,nu0,norm,c,zz,zcr,path,typo):
    
    pt_png=path+'//'+typo+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)
 
    lambda0=c/nu0
    
    if norm==1:
        zz=zz/zcr
        nu=nu/nu0
        nu_lim1=nu_lim1/nu0
        nu_lim2=nu_lim2/nu0
        xlabel='ν/ν0'
        ylabel='z/zcr'
    else:
        xlabel='ν, THz'
        ylabel='z, mm'
        nu=nu
        zz=zz
        
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes([0,0,1,1])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                    ['white','xkcd:bright orange','xkcd:chestnut'],
                                                    256)
    ax.set_title('η',fontsize=12)
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel(ylabel,fontsize=12)
    im=plt.imshow(eta,extent=[nu[0],nu[nu.shape[0]-1],zz[zz.shape[0]-1],zz[0]],aspect='auto',cmap = cmap, interpolation='none')
    ax.set_xlim(nu_lim1,nu_lim2)
    ax.locator_params(axis='x', nbins=5)
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax, shrink=2)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+str(1)+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf+'//'+str(1)+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()
    
    
def eta1d(eta1d,nu,nu_lim1,nu_lim2,nu0,norm,c,zz,zcr,path,typo):
    
    pt_png=path+'//'+typo+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    if norm==1:
        zz=zz/zcr
        nu=nu/nu0
        nu_lim1=nu_lim1/nu0
        nu_lim2=nu_lim2/nu0
        xlabel='z/zcr'
    else:
        xlabel='z, mm'
        nu=nu
        zz=zz

    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_subplot(111, label="1")
    plt.plot(zz, eta1d,'black')
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel('η',fontsize=12)
    # ax.set_xlim(t_lim1,t_lim2)
    ax.locator_params(axis='x', nbins=5)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+'E.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//'+'E.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def scatterMax2d(XX,YY,F,z,X,Nx,zcr,nu0,c,r_lim1,r_lim2,F_lim1,F_lim2,path,norm,typo):
    
    pt_png=path+'//'+typo+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    lambda0=c/nu0
    
    XX_mm=((XX-Nx/2)/Nx)*X
    YY_mm=((XX-Nx/2)/Nx)*X
    
    r=np.zeros(shape=XX.shape[0])
    for i in range(XX.shape[0]):
        r[i]=np.sqrt(XX_mm[i]**2+YY_mm[i]**2)
    
    if norm==1:
        r=r/lambda0
        r_lim1=r_lim1/lambda0
        r_lim2=r_lim2/lambda0
        xlabel='F#'
        ylabel='r/λ0'
        title = 'z/zcr= '+str(round(z/zcr,2))
    else:
        xlabel='F#'
        ylabel='z, mm'
        title = 'z= '+str(round(z,2))+' mm'
    

    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes([0,0,1,1])
    ax.set_title(title,fontsize=12)
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel(ylabel,fontsize=12)
    ax.scatter(F,r,color='xkcd:strong blue', zorder = 2, s = 4, marker = 'o')
    ax.set_xlim(F_lim1,F_lim2)
    ax.set_ylim(r_lim1,r_lim2)
    ax.locator_params(axis='x', nbins=5)
    # cax = fig55.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    # plt.colorbar(im, cax=cax, shrink=2)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.1)
    plt.savefig(pt_pdf+'//'+str(round(z,2))+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()
    

def scatterMax3d(XX,YY,F,z,X,Nx,zcr,nu0,c,x_lim1,x_lim2,y_lim1,y_lim2,z_lim1,z_lim2,path,norm,typo):
    
    pt_png=path+'//'+typo+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    lambda0=c/nu0
    
    XX_mm=((XX-Nx/2)/Nx)*X
    YY_mm=((XX-Nx/2)/Nx)*X
    
    if norm==1:
        XX_mm=XX_mm/lambda0
        YY_mm=YY_mm/lambda0
        x_lim1=x_lim1/lambda0
        x_lim2=x_lim2/lambda0
        y_lim1=y_lim1/lambda0
        y_lim2=y_lim2/lambda0
        xlabel='x/λ0'
        ylabel='y/λ0'
        title = 'z/zcr= '+str(round(z/zcr,2))
    else:
        XX_mm=XX_mm
        YY_mm=YY_mm
        F=F
        xlabel='x, mm'
        ylabel='y, mm'
        title = 'z= '+str(round(z,2))+' mm'
    
    # fig55 = plt.figure()
    fig = plt.figure(figsize=(2.5, 2.5))
    # ax = fig55.add_subplot(111, projection='3d')
    ax = Axes3D(fig)

    ax.set_title(title,fontsize=12)
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel(ylabel,fontsize=12)
    ax.set_zlabel('F#')

    ax.scatter(XX_mm,YY_mm,F,color='xkcd:strong blue', s = 4, marker = 'o')
    ax.set_xlim(x_lim1,x_lim2)
    ax.set_ylim(y_lim1,y_lim2)
    ax.set_zlim(z_lim1,z_lim2)
    
    ax.locator_params(axis='x', nbins=3)
    ax.locator_params(axis='y', nbins=3)
    # cax = fig55.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    # plt.colorbar(im, cax=cax, shrink=2)

    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+str(round(z,2))+'.png',bbox_inches = 'tight',pad_inches=0.4)
    plt.savefig(pt_pdf+'//'+str(round(z,2))+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def Ipart_z(Ipart,norm,c,zz,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//' + 'Ipart_z' + '//' + 'png'
    pt_pdf=path+'//'+typo+'//' + 'Ipart_z' + '//' + 'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    if norm==1:
        xlabel='z/zcr'
    else:
        xlabel='z, mm'
        zz=zz

    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_subplot(111, label="1")
    plt.plot(zz, Ipart,'black')
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel('I/Imax',fontsize=12)
    # ax.set_xlim(t_lim1,t_lim2)
    ax.locator_params(axis='x', nbins=5)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//' + str(prfx) + 'E.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//' + str(prfx) + 'E.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def FWHM_z(FWHM,norm,c,zz,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//' + 'FWHM' + '//' + 'png'
    pt_pdf=path+'//'+typo+'//' + 'FWHM' + '//' + 'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    if norm==1:
        xlabel='z/zcr'
    else:
        xlabel='z, mm'
        zz=zz

    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_subplot(111, label="1")
    plt.plot(zz, FWHM,'black')
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel('FWHM',fontsize=12)
    # ax.set_xlim(t_lim1,t_lim2)
    ax.locator_params(axis='x', nbins=5)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//' + str(prfx) + '_' + 'FWHM.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//' + str(prfx) + '_' + 'FWHM.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()  
    
    
def FWHM_z_6gr(FWHM,norm,c,zz,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//' + 'FWHM' + '//' + 'png'
    pt_pdf=path+'//'+typo+'//' + 'FWHM' + '//' + 'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    if norm==1:
        xlabel='z/zcr'
    else:
        xlabel='z, mm'
        zz=zz
    color_cycle = plt.rcParams['axes.prop_cycle']()
    
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_subplot(111, label="1")
    plt.plot(zz, FWHM[0,:],label='ρ=6.6', **next(color_cycle))
    plt.plot(zz, FWHM[1,:],label='ρ=10', **next(color_cycle))
    plt.plot(zz, FWHM[2,:],label='ρ=15', **next(color_cycle))
    plt.plot(zz, FWHM[3,:],label='ρ=20', **next(color_cycle))
    plt.plot(zz, FWHM[4,:],label='ρ=25', **next(color_cycle))
    plt.plot(zz, FWHM[5,:],label='ρ=30', **next(color_cycle))
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel('FWHM',fontsize=12)
    # ax.set_xlim(t_lim1,t_lim2)
    ax.set_ylim(0,30)
    ax.locator_params(axis='x', nbins=5)
    fig.set_dpi(100)
    
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

    fig.legend(lines, labels)
    
    plt.savefig(pt_png+'//' + str(prfx) + '_' + 'FWHM.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//' + str(prfx) + '_' + 'FWHM.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()
    

def ksi_z_6gr(ksis,norm,c,zz,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//' + 'ksis' + '//' + 'png'
    pt_pdf=path+'//'+typo+'//' + 'ksis' + '//' + 'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    if norm==1:
        xlabel='z/zcr'
    else:
        xlabel='z, mm'
        zz=zz
    color_cycle = plt.rcParams['axes.prop_cycle']()
    
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_subplot(111, label="1")
    plt.plot(zz, ksis[0,:],label='ρ=6.6', **next(color_cycle))
    plt.plot(zz, ksis[1,:],label='ρ=10', **next(color_cycle))
    plt.plot(zz, ksis[2,:],label='ρ=15', **next(color_cycle))
    plt.plot(zz, ksis[3,:],label='ρ=20', **next(color_cycle))
    plt.plot(zz, ksis[4,:],label='ρ=25', **next(color_cycle))
    plt.plot(zz, ksis[5,:],label='ρ=30', **next(color_cycle))
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel('ξ',fontsize=12)
    # ax.set_xlim(t_lim1,t_lim2)
    ax.set_ylim(0,1.1)
    ax.locator_params(axis='x', nbins=5)
    fig.set_dpi(100)
    
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

    fig.legend(lines, labels)
    
    plt.savefig(pt_png+'//' + str(prfx) + '_' + 'FWHM.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//' + str(prfx) + '_' + 'FWHM.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()


def graph1d(ar,path,typo,prfx):
    
    pt_png=path+'//'+typo+'//' + 'En' + '//' + 'png'
    pt_pdf=path+'//'+typo+'//' + 'En' + '//' + 'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    xlabel='iter'

    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_subplot(111, label="1")
    plt.plot(ar,'black')
    plt.xlabel(xlabel,fontsize=12)
    plt.ylabel('FWHM',fontsize=12)
    # ax.set_xlim(t_lim1,t_lim2)
    ax.locator_params(axis='x', nbins=5)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//' + str(prfx) + '_' + 'En.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//' + str(prfx) + '_' + 'En.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()
    

def traj(traj,X,Y,Nx,Ny,z,typo,path,prfx):
    
    pt_png=path+'//'+typo+'//'+'traj'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'traj'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
          
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes([0,0,1,1])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                    ['xkcd:strong blue','white','xkcd:red'],
                                                    256)
    ax.set_title('trajectories',fontsize=12)
    plt.xlabel('x, mm',fontsize=12)
    plt.ylabel('y, mm',fontsize=12)
    im=plt.imshow(traj,extent=[x[0],x[Nx-1],y[Ny-1],y[0]],aspect='auto',cmap = cmap, interpolation='none')
    plt.clim(-1, 1)
    ax.locator_params(axis='x', nbins=5)
    ax.locator_params(axis='y', nbins=5)
    cax = fig.add_axes([ax.get_position().x1+0.02,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax, shrink=2)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+ str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//'+ str(prfx) + '_' + str(round(z,2))+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()
    
    
def traj_video(traj,nu,X,Y,Nx,Ny,z,typo,path,prfx):
    
    pt_png=path+'//'+typo+'//'+'traj_video'+'//'+ str(round(z,2))+'//'+'png'
    pt_pdf=path+'//'+typo+'//'+'traj_video'+'//'+ str(round(z,2))+'//'+'pdf'
    
    if not os.path.isdir(pt_png):
        os.makedirs(pt_png)
        
    if not os.path.isdir(pt_pdf):
        os.makedirs(pt_pdf)

    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
          
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes([0,0,1,1])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                    ['xkcd:strong blue','white','xkcd:red'],
                                                    256)
    # ax.set_title('trajectories',fontsize=12)
    ax.set_title('nu= ' + str(round((nu),2)),fontsize=12)
    plt.xlabel('x, mm',fontsize=12)
    plt.ylabel('y, mm',fontsize=12)
    im=plt.imshow(traj,extent=[x[0],x[Nx-1],y[Ny-1],y[0]],aspect='auto',cmap = cmap, interpolation='none')
    plt.clim(-1, 1)
    ax.locator_params(axis='x', nbins=5)
    ax.locator_params(axis='y', nbins=5)
    cax = fig.add_axes([ax.get_position().x1+0.02,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax, shrink=2)
    fig.set_dpi(80)
    plt.savefig(pt_png+'//'+ str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight')
    plt.savefig(pt_pdf+'//'+ str(prfx) + '_' + str(round(z,2))+'.pdf',bbox_inches = 'tight')
    plt.cla() 
    plt.clf() 
    plt.close('all')   
    plt.close(fig)
    gc.collect()

    
def charges(charges,nu,X,Y,Nx,Ny,z,typo,path,prfx):
   
   pt_png=path+'//'+typo+'//'+'charges'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'png'
   pt_pdf=path+'//'+typo+'//'+'charges'+'//'+ str(prfx) + '_' + str(round(z,2))+'//'+'pdf'
   
   if not os.path.isdir(pt_png):
       os.makedirs(pt_png)
       
   if not os.path.isdir(pt_pdf):
       os.makedirs(pt_pdf)
         
   fig = plt.figure(figsize=(2, 2))
   ax = fig.add_axes([0,0,1,1])
   cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                   ['xkcd:strong blue','white','xkcd:red'],
                                                   256)
   ax.set_title('charges',fontsize=12)
   plt.xlabel('ν, THz',fontsize=12)
   im=plt.imshow(charges,extent=[nu[0],nu[nu.shape[0]-1],0,29],aspect=0.005,cmap = cmap, interpolation='none')
   plt.clim(-1, 1)
   ax.locator_params(axis='x', nbins=10)
   ax.locator_params(axis='y', nbins=5)
   # plt.rc('xtick',labelsize=8)
   ax.set_yticks([])
   # cax = fig.add_axes([ax.get_position().x1+0.02,ax.get_position().y0,0.02,ax.get_position().height])
   # plt.colorbar(im, cax=cax, shrink=2)
   fig.set_dpi(150)
   plt.savefig(pt_png+'//'+ str(prfx) + '_' + str(round(z,2))+'.png',bbox_inches = 'tight')
   plt.savefig(pt_pdf+'//'+ str(prfx) + '_' + str(round(z,2))+'.pdf',bbox_inches = 'tight')
   plt.cla() 
   plt.clf() 
   plt.close('all')   
   plt.close(fig)
   gc.collect()   
   
    
    
    
    
    
    
    