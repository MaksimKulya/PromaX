import matplotlib.pyplot as plt
import math
import cmath
import pylab
from matplotlib import mlab
from PIL import Image
import numpy as np
import scipy.constants
import numba as nb

import AS
import C

pi=math.pi


def Prop(GG0,nu,X,Y,nmedia,c,z,nu_cr):
    
    nu_cr_ind = np.argmin(np.abs(nu - nu_cr))
    
    if z==0:
        GG=GG0
    else:
        GG0_C=GG0[0:nu_cr_ind,:,:]
        GG_C=C.C(GG0_C,nu[0:nu_cr_ind],X,Y,nmedia[0:nu_cr_ind],c,z)

        GG0_AS=GG0[nu_cr_ind:GG0.shape[0],:,:]
        GG_AS=AS.AS(GG0_AS,nu[nu_cr_ind:GG0.shape[0]],X,Y,nmedia[nu_cr_ind:GG0.shape[0]],c,z)
    
        GG=np.concatenate([GG_C, GG_AS], 0)
    
    return GG



def Propz(GG0,nu,X,Y,nmedia,c,z,nu_cr):
    
    nu_cr_ind = np.argmin(np.abs(nu - nu_cr))
    
    if z==0:
        GGz=AS.ASz(GG0,nu,X,Y,nmedia,c,z)
    else:
        GG0_Cz=GG0[0:nu_cr_ind,:,:]
        GG_Cz=C.Cz(GG0_Cz,nu[0:nu_cr_ind],X,Y,nmedia[0:nu_cr_ind],c,z)

        GG0_ASz=GG0[nu_cr_ind:GG0.shape[0],:,:]
        GG_ASz=AS.ASz(GG0_ASz,nu[nu_cr_ind:GG0.shape[0]],X,Y,nmedia[nu_cr_ind:GG0.shape[0]],c,z)
    
        GGz=np.concatenate([GG_Cz, GG_ASz], 0)
    
    return GGz


def Propz_polar(GG0x,GG0y,nu,X,Y,nmedia,c,z,nu_cr):
    
    nu_cr_ind = np.argmin(np.abs(nu - nu_cr))
    
    if z==0:
        GGz=AS.ASz_polar(GG0x,GG0y,nu,X,Y,nmedia,c,z)
    else:
        GG0x_Cz=GG0x[0:nu_cr_ind,:,:]
        GG0y_Cz=GG0y[0:nu_cr_ind,:,:]
        GG_Cz=C.Cz_polar(GG0x_Cz,GG0y_Cz,nu[0:nu_cr_ind],X,Y,nmedia[0:nu_cr_ind],c,z)

        GG0x_ASz=GG0x[nu_cr_ind:GG0x.shape[0],:,:]
        GG0y_ASz=GG0y[nu_cr_ind:GG0x.shape[0],:,:]
        GG_ASz=AS.ASz_polar(GG0x_ASz,GG0y_ASz,nu[nu_cr_ind:GG0x.shape[0]],X,Y,nmedia[nu_cr_ind:GG0x.shape[0]],c,z)
    
        GGz=np.concatenate([GG_Cz, GG_ASz], 0)
    
    return GGz


