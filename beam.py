import matplotlib.pyplot as plt
import math
import cmath
import pylab
from matplotlib import mlab
from PIL import Image
import numpy as np
import scipy.constants
import numba as nb


def beam(G,Nx,Ny):

    AM = np.zeros(shape=(G.shape[0], Nx, Ny))
    PH = np.zeros(shape=(G.shape[0], Nx, Ny))
    G_beam = np.zeros(shape=(G.shape[0], Nx, Ny), dtype=complex)

    am_ones = np.ones(shape=(Nx, Ny))

    for k in range(G.shape[0]):
        AM[k, :, :] = am_ones * np.abs(G[k])
        PH[k, :, :] = np.angle(G[k])
        G_beam[k, :, :] = AM[k, :, :] * np.exp(1j * PH[k, :, :])

    return G_beam


