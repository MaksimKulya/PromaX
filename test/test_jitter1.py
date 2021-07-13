import matplotlib
import matplotlib.pyplot as plt
import numpy as np






fig1, ax1 = plt.subplots()
c = ax1.imshow(np.random.random((10, 10)))
cbar_ax = fig1.add_axes([0.1, 0.1, 0.05, 0.8])
# new ax with dimensions of the colorbar

cbar = fig1.colorbar(c, cax=cbar_ax)
plt.savefig('D:\\SCIENCE\\2021\\RNF\\1\\1.png',bbox_inches = 'tight',dpi=200/fig1.get_size_inches()[1]) 




fig2, ax2 = plt.subplots()
c = ax2.imshow(-1*np.random.random((10, 10)))
cbar_ax = fig2.add_axes([0.1, 0.1, 0.05, 0.8])
# new ax with dimensions of the colorbar

cbar = fig2.colorbar(c, cax=cbar_ax)
plt.savefig('D:\\SCIENCE\\2021\\RNF\\1\\2.png',bbox_inches = 'tight',dpi=200/fig2.get_size_inches()[1])


plt.show()
