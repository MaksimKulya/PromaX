import matplotlib.pyplot as plt
import math
import cmath
import pylab
from matplotlib import mlab
from PIL import Image
import numpy as np
import scipy.constants
import time

# read amplitude and phase masks
def read_mask(pathAM,pathPH):
    im1 = Image.open(pathAM)
    AM = np.array(im1)
    AM = AM/np.amax(AM)

    im2 = Image.open(pathPH)
    PH = np.array(im2)
    PH = PH/np.amax(PH)
    
    return AM,PH

