import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from itertools import permutations 
import math

path_qwp = 'D:\\SCIENCE\\2021\\RNF\\Mar\\plates\\waveplate_quarter_05-13.txt'

wp = np.loadtxt(path_qwp, comments="#", delimiter=",", unpack=False)

h=wp[:,0]

angles = wp[:,1]



perm1 = permutations(h)
perm2 = permutations(h)

hh = np.zeros(shape=math.factorial(len(h)), dtype=float)  


perm1=np.array(list(perm1))


for i in range(perm1.shape[0]):
    h = perm1[i,:]