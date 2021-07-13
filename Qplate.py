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
import numba as nb

import Jones
import Qplate_core

pi=math.pi
c = scipy.constants.c * 1e-9 #speed of light in mm/ps

# 8 sectors 
def Qplate_sr(path,GG0,X,Y,Nx,Ny,nu,GG0x,GG0y,del_b,Lf):
    # 7 q/2
    wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
    h=wp[:,0]
    h=h[::-1]
    
    angles = wp[:,1]
    angles = angles * pi / 180
    angles=angles[::-1]
    
    n_o=2.115
    n_e=2.165
    
    del_b = del_b * pi / 180
    # -----------------------------------------------------------------------------
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
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
    # -----------------------------------------------------------------------------
    GG0x_mod15 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod15 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
    GG0x_mod9 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod9 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
    GG0x_mod4 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod4 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
    GG0x_mod8 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod8 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
    GG0x_mod3 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod3 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
    GG0x_mod7 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod7 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
    GG0x_mod2 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod2 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    
    GG0x_mod6 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)  
    GG0y_mod6 = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
      
    # segment 1+5 
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q]+0+del_b)
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
        GG0x_mod15[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*phi15*rho2
        GG0y_mod15[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*phi15*rho2
    
    # segment 9        
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q]-157.5*pi/180+del_b)
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
        GG0x_mod9[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*phi9*rho2
        GG0y_mod9[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*phi9*rho2
    
    # segment 4  
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q]+45*pi/180+del_b)
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
        GG0x_mod4[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*phi4*rho2
        GG0y_mod4[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*phi4*rho2
    
    # segment 8  
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q]-112.5*pi/180+del_b)
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
        GG0x_mod8[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*phi8*rho2
        GG0y_mod8[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*phi8*rho2
        
    # segment 3  
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q]+90*pi/180+del_b)
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
        GG0x_mod3[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*phi3*rho2
        GG0y_mod3[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*phi3*rho2
                
    # segment 7  
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q]-67.5*pi/180+del_b)
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
        GG0x_mod7[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*phi7*rho2
        GG0y_mod7[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*phi7*rho2
               
    # segment 2
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q]+135*pi/180+del_b)
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
        GG0x_mod2[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*phi2*rho2
        GG0y_mod2[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*phi2*rho2
              
    # segment 6
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q]-22.5*pi/180+del_b)
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
        GG0x_mod6[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*phi6*rho2
        GG0y_mod6[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*phi6*rho2
    
    GG0x_qplate = GG0x_mod15 + GG0x_mod9 + GG0x_mod4 + GG0x_mod8 + GG0x_mod3 + GG0x_mod7 + GG0x_mod2 + GG0x_mod6
    GG0y_qplate = GG0y_mod15 + GG0y_mod9 + GG0y_mod4 + GG0y_mod8 + GG0y_mod3 + GG0y_mod7 + GG0y_mod2 + GG0y_mod6
    
    return GG0x_qplate, GG0y_qplate


def Qplate_360r(path,GG0,X,Y,Nx,Ny,nu,GG0x,GG0y,del_b,Lf):
    # 7 q/2
    wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
    h=wp[:,0]
    h=h[::-1]
    
    angles = wp[:,1]
    angles = angles * pi / 180
    angles=angles[::-1]
    
    n_o=2.115
    n_e=2.165
    
    del_b = del_b * pi / 180
    # -----------------------------------------------------------------------------
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
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
        PHI[k,:,:]=(np.max(phix[k,:,:])-np.min(phix[k,:,:]))/Lf   
        PHI[k,:,:] = PHI[k,:,:] * sectorx[k,:,:]
    
    # take 1st non zero element and put it into new 1d array    
    PHI1 = np.zeros(shape=(N3d), dtype=float)
    for k in range(N3d):    
        for i in range(Nx):
                for j in range(Ny):    
                    if PHI[k,i,j] != 0:
                        PHI1[k] = PHI[k,i,j]
                        break
                    
# -----------------------------------------------------------------------------    

    GG0x_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    GG0y_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex) 
    
    GG0x_qplate = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    GG0y_qplate = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
        
    for w in range(N3d):
        
        for k in range(len(nu)):    
            jonesmatrix = np.eye(2, dtype=complex)    
            for q in range(len(h)):     
            
                delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
                jonesmatrice = Jones.get_jonesmatrix(delta, angles[q] + PHI1[w] + del_b)
                jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
                
            GG0x_mod[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*sectorx[w,:,:]*rho2
            GG0y_mod[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*sectorx[w,:,:]*rho2
                   
        GG0x_qplate = GG0x_qplate + GG0x_mod 
        GG0y_qplate = GG0y_qplate + GG0y_mod 
    
    return GG0x_qplate, GG0y_qplate


def Qplate_360rnb(path,GG0,X,Y,Nx,Ny,nu,GG0x,GG0y,del_b,Lf):
    # 7 q/2
    wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
    h=wp[:,0]
    h=h[::-1]
    
    angles = wp[:,1]
    angles = angles * pi / 180
    angles=angles[::-1]
    
    n_o=2.115
    n_e=2.165
    
    del_b = del_b * pi / 180
    # -----------------------------------------------------------------------------
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
    N3d=100

    (PHI1,sectorx,rho2) = Qplate_core.loopa(N3d,Nx,Ny,x,y,X,Lf)
        
    GG0x_qplate,GG0y_qplate = Qplate_core.loopa2(GG0,N3d,nu,h,n_e,n_o,angles,PHI1,del_b,GG0x,GG0y,sectorx,rho2)
    
    return GG0x_qplate, GG0y_qplate


def Qplate_360_1layer(path,GG0,X,Y,Nx,Ny,nu,GG0x,GG0y,del_b,h,angles):
    # 7 q/2

    # h=h[::-1]
    # angles=angles[::-1]
    
    n_o=2.115
    n_e=2.165
    
    del_b = del_b * pi / 180
    # -----------------------------------------------------------------------------
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
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
    
    # take 1st non zero element and put it into new 1d array    
    PHI1 = np.zeros(shape=(N3d), dtype=float)
    for k in range(N3d):    
        for i in range(Nx):
                for j in range(Ny):    
                    if PHI[k,i,j] != 0:
                        PHI1[k] = PHI[k,i,j]
                        break
# -----------------------------------------------------------------------------    

    GG0x_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    GG0y_mod = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex) 
    
    GG0x_qplate = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
    GG0y_qplate = np.zeros(shape=(GG0.shape[0], GG0.shape[1], GG0.shape[2]), dtype=complex)
        
    for w in range(N3d):
        
        for k in range(len(nu)):    
            jonesmatrix = np.eye(2, dtype=complex)    
            for q in range(len(h)):     
            
                delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
                jonesmatrice = Jones.get_jonesmatrix(delta, angles[q] + PHI1[w] + del_b)
                jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
                
            GG0x_mod[k,:,:] = (jonesmatrix[0,0]*GG0x[k,:,:] + jonesmatrix[0,1]*GG0y[k,:,:])*sectorx[w,:,:]*rho2
            GG0y_mod[k,:,:] = (jonesmatrix[1,0]*GG0x[k,:,:] + jonesmatrix[1,1]*GG0y[k,:,:])*sectorx[w,:,:]*rho2
                   
        GG0x_qplate = GG0x_qplate + GG0x_mod 
        GG0y_qplate = GG0y_qplate + GG0y_mod 
    
    return GG0x_qplate, GG0y_qplate


def efficency(path,GG0,X,Y,Nx,Ny,nu):
    # 7 q/2
    wp = np.loadtxt(path, comments="#", delimiter=",", unpack=False)
    
    h=wp[:,0]
    
    angles = wp[:,1]
    angles = angles * pi / 180
    
    n_o=2.115
    n_e=2.165

    retard = np.zeros(shape=(len(nu),1))
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q])
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
        retard[k] = Jones.get_delay(jonesmatrix)/pi

    
    return retard


def efficency_1layer(path,GG0,X,Y,Nx,Ny,nu,h,angles):
    # 7 q/2
    
    n_o=2.115
    n_e=2.165

    retard = np.zeros(shape=(len(nu),1))
    for k in range(len(nu)):    
        jonesmatrix = np.eye(2, dtype=complex)    
        for q in range(len(h)):     
        
            delta = 2*pi*h[q]*(n_e-n_o)*nu[k]/c
            jonesmatrice = Jones.get_jonesmatrix(delta, angles[q])
            jonesmatrix=np.matmul(jonesmatrix,jonesmatrice)
            
        retard[k] = Jones.get_delay(jonesmatrix)/pi

    
    return retard








