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

N=16

nucrop=4

a=np.array([[1,1],[1,1]])

aa=np.zeros(shape=(N,2,2))

for k in range(N):
    aa[k,:,:]=a*k
    
aa_crop=aa[(N//2)+1:(N//2)+1+nucrop,:,:]    

GG_ansig = np.pad(aa_crop, ((N//2+1, N//2-aa_crop.shape[0]-1), (0, 0), (0, 0)), mode='constant')

# aa_cr = aa_zpd[:,aa.shape[1]//2:aa.shape[1]//2+aa.shape[1],aa.shape[1]//2:aa.shape[1]//2+aa.shape[1]]




