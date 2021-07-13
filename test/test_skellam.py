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
from scipy.special import *

from scipy.stats import skellam
from scipy.stats import norm


xx = np.linspace(-10,10,501)
# -----------------------------------------------------------------------------
# Modified Bessel 1st order
J=((1j**(-abs(xx)))*jv(abs(xx),1j*1)).real  

# fig1, ax = plt.subplots(1, 1)
# ax.plot(xx, J, 'bo', ms=2, label='skellam pmf')
# -----------------------------------------------------------------------------
# Skellam by eq.
I0=1
Gamma=0.9
mu1=I0*(1+Gamma)/2
mu2=I0*(1-Gamma)/2

skellam_fx=np.exp(-mu1-mu2)*((mu1/mu2)**(xx/2))*J*2*np.sqrt(mu1*mu2)

mu_skellam=mu1-mu2
sigma_skellam=np.sqrt(mu1+mu2)


# -----------------------------------------------------------------------------
# Gauss
mu=mu1-mu2
sigma=np.sqrt(mu1+mu2)

gauss_fx = norm(mu, sigma)
# -----------------------------------------------------------------------------

fig, ax = plt.subplots(1, 1)
ax.plot(xx, skellam_fx, 'ro', ms=1, label='skellam pmf')
ax.plot(xx, gauss_fx.pdf(xx),'bo',ms=1)


ax.text(-7, 0.3, "Gauss (blue):" + "\n"
        + "mu= " + str(np.round(mu,2)) + "\n" 
        + "sigma= " + str(np.round(sigma,2)) + "\n", ha="center", va="center", size=12)

ax.text(7, 0.3, "Skellam(red):" + "\n" +
        "I0= " + str(I0) + "\n" +
        "Gamma= " + str(Gamma) + "\n" +
        "mu1= " + str(np.round(mu1,2)) + "\n"+
        "mu2= " + str(np.round(mu2,2)) + "\n"+ 
        "mu= " + str(np.round(mu_skellam,2)) + "\n"+
        "sigma= " + str(np.round(sigma_skellam,2)), ha="center", va="center", size=12)






# fig2, ax = plt.subplots(1, 1)
# ax.plot(xx, skellam.pmf(xx, mu1, mu2), 'go', ms=2, label='skellam pmf')
# # ax.vlines(xx, 0, skellam.pmf(xx, mu1, mu2), colors='b', lw=5, alpha=0.5)
# # ax.plot(xx, rv.pdf(xx))

# ax.text(-7, 0.35, "Gauss:" + "\n"
#         + "mu= " + str(mu) + "\n" 
#         + "sigma= " + str(sigma) + "\n", ha="center", va="center", size=12)

# ax.text(7, 0.35, "Skellam:" + "\n" +
#         "I0= " + str(I0) + "\n" 
#         + "Gamma= " + str(Gamma) + "\n"
#         + "sigma= " + str(round(sigma_skellam,2)), ha="center", va="center", size=12)













