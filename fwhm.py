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
from scipy.interpolate import interp1d

def peak(x, c):
    return np.exp(-np.power(x - c, 2) / 16.0)

def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def half_max_x(x, y, level):
    half = max(y)/level
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    hmx = [lin_interp(x, y, zero_crossings_i[0], half),
            lin_interp(x, y, zero_crossings_i[1], half)]
    fwhm = hmx[1] - hmx[0]
    return fwhm



# def fwhm(x, y):
    
#     f = interp1d(x, y, kind='cubic')
#     x = np.linspace(x[0], x[-1], len(x)*5, endpoint=True)
#     y = f(x)

#     max_y = max(y)  # Find the maximum y value
#     ii = [i for i in range(len(y)) if y[i] > max_y/2]

#     fwhm = x[max(ii)]-x[min(ii)]

#     return fwhm


def get_fwhm(x, GG, y0):
    
    GG2d=abs(GG).sum(axis=0)
    GG1d = GG2d[:,y0]
    
    f = interp1d(x, GG1d, kind='cubic')
    x = np.linspace(x[0], x[-1], len(x)*5, endpoint=True)
    GG1d = f(x)

    max_GG1d = max(GG1d)  # Find the maximum y value
    ii = [i for i in range(len(GG1d)) if GG1d[i] > max_GG1d/2]

    fwhm = x[max(ii)]-x[min(ii)]

    return fwhm


# def get_fwhm2(x, GG, GG0, y0):
    
#     GG2d=abs(GG).sum(axis=0)
#     GG1d = GG2d[:,y0]
    
#     f = interp1d(x, GG1d, kind='cubic')
#     x = np.linspace(x[0], x[-1], len(x)*5, endpoint=True)
#     GG1d = f(x)

#     max_GG1d = max(GG1d)  # Find the maximum y value
#     ii = [i for i in range(len(GG1d)) if GG1d[i] > max_GG1d/2]

#     fwhm = x[max(ii)]-x[min(ii)]

#     GG02d=abs(GG0).sum(axis=0)
           
#     a = GG2d[min(ii)//5:max(ii)//5,min(ii)//5:max(ii)//5]
#     a = a.sum()
    
#     b = GG02d.sum()
    
#     ksi = a/b

#     return fwhm, ksi


def En_ondet(X, Nx, size_det, GG, GG0):
    
    det = int(size_det/X*Nx)
    
    GG_ondet = GG[:,Nx//2-det//2:Nx//2+det//2,Nx//2-det//2:Nx//2+det//2]  
   
    ksi = abs(GG_ondet)**2
    init = abs(GG0)**2

    ksi = ksi.sum()
    init = init.sum()
    
    ksi = ksi/init

    return ksi
















