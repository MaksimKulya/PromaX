import matplotlib.pyplot as plt
import math
import cmath
import pylab
from matplotlib import mlab
from PIL import Image
import numpy as np
import scipy.constants

pi=math.pi

def SLM(pathSLM,nu,index_800nm,Nx,Ny):

    slm = Image.open(pathSLM)
    SLM = np.array(slm)
    SLM = SLM/np.amax(SLM)    
    phi=(SLM*2-1)*pi*(-1)

    phi_SLM=np.zeros(shape=(nu.shape[0], Nx, Ny), dtype=complex)
    G_SLM=np.zeros(shape=(nu.shape[0], Nx, Ny), dtype=complex)
    
    for k in range(nu.shape[0]):
        phi_SLM[k,:,:]=(phi*nu[k]/nu[index_800nm])
        G_SLM[k,:,:]=np.exp(1j*phi_SLM[k,:,:])

    return G_SLM


