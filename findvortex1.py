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

import numba as nb

import gc
import os


from numpy import unravel_index


# numin=20
# numax=300
# deepsearch=1000
# cropxy=22
# sizePs=10

@nb.njit
def fnd(GG,Nx,Ny,X,Y,nu,numin,numax,deepsearch,cropxy,sizePs,z,path,typo,prfx):

    traj = np.zeros(shape=(Nx,Ny), dtype=nb.types.float64)
    
    nu=nu[numin:numax]
    charges = np.zeros(shape=(30,len(nu)), dtype=nb.types.float64)
    
    for k in range(numin,numax):
    
        GG1 = GG[k,:,:]
        
        pi=math.pi
        
        # find singular by ampl
        
        A = np.power(np.absolute(GG1), 1/15)
        
        ch=0
        count=0

        while ch != 1 and count<deepsearch:
            result = np.where(A == np.amin(A))
            ind = list(zip(result[0], result[1]))
            
            ind = np.asarray(ind, dtype=nb.types.float64)
            indx=int(ind[0,0])
            indy=int(ind[0,1])
            

            if indx<cropxy or indx>Nx-cropxy or indy<cropxy or indy>Ny-cropxy:
                A[indx,indy]=100
                count=count+1
                continue
                
            P = np.angle(GG1)
            
            Ps = P[indx-sizePs:indx+sizePs,indy-sizePs:indy+sizePs]
            
            # find phase rotation
            
            NNx = 2*sizePs
            NNy = 2*sizePs
            
            XX = X * (NNx/Nx)
            YY = Y * (NNy/Ny)
            
            phi1 = np.zeros(shape=(NNx,NNy), dtype=nb.types.float64)
            
            N3d=1000
            
            ph=np.linspace(-pi,pi,N3d)
            r=XX/5
            
            xcoor=np.zeros(shape=(N3d), dtype=nb.types.float64)
            ycoor=np.zeros(shape=(N3d), dtype=nb.types.float64)
            
            for i,phi in enumerate(ph):    
                xcoor[i] = int(r * np.cos(phi) * NNx/XX)
                ycoor[i] = int(r * np.sin(phi) * NNx/XX)
            
            PPs=np.zeros(shape=(N3d), dtype=nb.types.float64)
            
            for i in range(len(xcoor)):    
                PPs[i] = Ps[int(xcoor[i])+len(Ps)//2, int(ycoor[i])+len(Ps)//2]
            
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
            else:A[indx,indy]=100
            
            count=count+1
            
            
    # visual.traj(traj,X,Y,Nx,Ny,z,typo,path,prfx)
    # visual.charges(charges,nu,X,Y,Nx,Ny,z,typo,path,prfx)





