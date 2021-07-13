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
from itertools import permutations 

perm2 = np.zeros(shape=(100, 100), dtype=float) 

np.set_printoptions(threshold=np.inf, linewidth=np.inf)  # turn off summarization, line-wrapping
with open("perm2.txt", 'w') as f:
    f.write(np.array2string(perm2, separator=', '))
