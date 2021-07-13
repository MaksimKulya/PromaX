import matplotlib
import matplotlib.pyplot as plt
import math
import cmath
from PIL import Image
import numpy as np
import scipy.constants
import time


N=5
x=np.linspace(1,N,N)
y=np.linspace(1,N,N)

nu=np.linspace(0,100,11)

# Xm,NUm,Ym = np.meshgrid(x,nu,y)
Xm,Ym = np.meshgrid(x,y)





