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



a=np.array([[1,1],[1,1]])

aa=np.zeros(shape=(8,2,2))

for k in range(4):
    aa[k,:,:]=a*(k+1)

aa_zpd = np.pad(aa, ((0, 0), (aa.shape[1]//2, aa.shape[1]//2), (aa.shape[1]//2, aa.shape[1]//2)), mode='constant')

aa_cr = aa_zpd[:,aa.shape[1]//2:aa.shape[1]//2+aa.shape[1],aa.shape[1]//2:aa.shape[1]//2+aa.shape[1]]




