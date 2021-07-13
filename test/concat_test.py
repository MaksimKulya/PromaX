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


Nt=7
Nx=2
Ny=2

nu_cr_ind=3

a=np.array([[1,1],[1,1]])

aa=np.zeros(shape=(Nt,Nx,Ny))

for k in range(Nt):
    aa[k,:,:]=a*(k+1)

a1=aa[0:nu_cr_ind,:,:]

a2=aa[nu_cr_ind:aa.shape[0],:,:]


b=np.concatenate([a1, a2], 0)


x=[1,2,3,4]

xx=x.index(3)

print(aa.shape[0])