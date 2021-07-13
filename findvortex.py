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

import visual

import gc
import os


from numpy import unravel_index


# numin=20
# numax=300
# deepsearch=1000
# cropxy=22
# sizePs=10

def fnd(GG,Nx,Ny,X,Y,nu,numin,numax,deepsearch,cropxy,sizePs,z,path,typo,prfx):
    
    pi=math.pi

    traj = np.zeros(shape=(Nx,Ny), dtype=float)
    
    nu=nu[numin:numax]
    charges = np.zeros(shape=(30,len(nu)), dtype=float)
    
    for k in range(numin,numax):
    
        GG1 = GG[k,:,:]
        
        # find singular by ampl
        
        A = np.power(abs(GG1), 1/15)
        P = np.angle(GG1)
        
        NNx = 2*sizePs
        XX = X * (NNx/Nx)
        N3d=64
        ph=np.linspace(-pi,pi,N3d)
        r=XX/2
        
        ch=0
        count=0

        while ch != 1 and count<deepsearch:
            result = np.where(A == np.amin(A))
            ind = list(zip(result[0], result[1]))
            
            ind = np.asarray(ind, dtype=np.float)
            indx=int(ind[0,0])
            indy=int(ind[0,1])
            
            if indx<cropxy or indx>Nx-cropxy or indy<cropxy or indy>Ny-cropxy:
                A[indx,indy]=1000
                count=count+1
                continue
            
            Ps = P[indx-sizePs:indx+sizePs,indy-sizePs:indy+sizePs]
            
            # find phase rotation         
            PPs=np.zeros(shape=(N3d), dtype=float)
            
            for i,phi in enumerate(ph):    
                PPs[i] = Ps[int(r*np.cos(phi)*NNx/XX)+len(Ps)//2, int(r*np.sin(phi)*NNx/XX)+len(Ps)//2]
            
            PPs=np.unwrap(PPs)
            
            charge = (PPs[-1]-PPs[0])/2/pi

            # draw trajectory
            
            if (charge==1 and PPs[0]<PPs[N3d//4]<PPs[N3d//2]<PPs[N3d//2+N3d//4]<PPs[N3d-1]):  
                traj[indx,indy]=1
                charges[:,k-numin]=1
                ch=1
            elif (charge==-1 and PPs[0]>PPs[N3d//4]>PPs[N3d//2]>PPs[N3d//2+N3d//4]>PPs[N3d-1]):
                traj[indx,indy]=-1
                charges[:,k-numin]=-1
                ch=1
            else:A[indx,indy]=1000
            
            count=count+1
            
            
    visual.traj(traj,X,Y,Nx,Ny,z,typo,path,prfx)
    visual.charges(charges,nu,X,Y,Nx,Ny,z,typo,path,prfx)


def fnd_video(GG,Nx,Ny,X,Y,nu,numin,numax,deepsearch,cropxy,sizePs,z,path,typo,prfx):
    
    pi=math.pi

    traj = np.zeros(shape=(Nx,Ny), dtype=float)
    
    nu=nu[numin:numax]
    charges = np.zeros(shape=(30,len(nu)), dtype=float)
    
    for k in range(numin,numax):
    
        GG1 = GG[k,:,:]
        
        # find singular by ampl
        
        A = np.power(abs(GG1), 1/15)
        P = np.angle(GG1)
        
        NNx = 2*sizePs
        XX = X * (NNx/Nx)
        N3d=64
        ph=np.linspace(-pi,pi,N3d)
        r=XX/2
        
        ch=0
        count=0

        while ch != 1 and count<deepsearch:
            result = np.where(A == np.amin(A))
            ind = list(zip(result[0], result[1]))
            
            ind = np.asarray(ind, dtype=np.float)
            indx=int(ind[0,0])
            indy=int(ind[0,1])
            
            if indx<cropxy or indx>Nx-cropxy or indy<cropxy or indy>Ny-cropxy:
                A[indx,indy]=1000
                count=count+1
                continue
            
            Ps = P[indx-sizePs:indx+sizePs,indy-sizePs:indy+sizePs]
            
            # find phase rotation         
            PPs=np.zeros(shape=(N3d), dtype=float)
            
            for i,phi in enumerate(ph):    
                PPs[i] = Ps[int(r*np.cos(phi)*NNx/XX)+len(Ps)//2, int(r*np.sin(phi)*NNx/XX)+len(Ps)//2]
            
            PPs=np.unwrap(PPs)
            
            charge = (PPs[-1]-PPs[0])/2/pi

            # draw trajectory
            
            if (charge==1 and PPs[0]<PPs[N3d//4]<PPs[N3d//2]<PPs[N3d//2+N3d//4]<PPs[N3d-1]):  
                traj[indx,indy]=1
                charges[:,k-numin]=1
                visual.traj_video(traj,nu[k],X,Y,Nx,Ny,z,typo,path,k)
                ch=1
            elif (charge==-1 and PPs[0]>PPs[N3d//4]>PPs[N3d//2]>PPs[N3d//2+N3d//4]>PPs[N3d-1]):
                traj[indx,indy]=-1
                charges[:,k-numin]=-1
                visual.traj_video(traj,nu[k],X,Y,Nx,Ny,z,typo,path,k)
                ch=1
            else:A[indx,indy]=1000
            
            count=count+1

    visual.traj(traj,X,Y,Nx,Ny,z,typo,path,prfx)
    visual.charges(charges,nu,X,Y,Nx,Ny,z,typo,path,prfx)
