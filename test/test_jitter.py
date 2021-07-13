import matplotlib
import matplotlib.pyplot as plt
import numpy as np

array1=np.ones(shape=(50,50))
array2=-np.ones(shape=(50,50))
# -----------------------------------------------------------------------------
fig1 = plt.figure(figsize=(3, 3))
fig1.set_size_inches(2,2)
ax1=fig1.add_subplot(111, label="1")
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                ['white','xkcd:bright orange','xkcd:chestnut'],
                                                256)
plt.xlabel('x',fontsize=12)
plt.ylabel('y',fontsize=12)

im=plt.imshow(array1,aspect='equal',cmap = cmap)
cax1 = fig1.add_axes([1, 0.1, 0.05, 0.8])
# cax1 = fig1.add_axes([ax1.get_position().x1+0.01,ax1.get_position().y0,0.02,ax1.get_position().height])

fig1.colorbar(im, cax=cax1, shrink=2)
plt.savefig('D:\\SCIENCE\\2021\\RNF\\1\\1.png',bbox_inches = 'tight',dpi=200/fig1.get_size_inches()[1])           
            
# -----------------------------------------------------------------------------
fig2 = plt.figure(figsize=(3, 3))
fig2.set_size_inches(2,2)
ax2=fig2.add_subplot(111, label="1")
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                ['white','xkcd:bright orange','xkcd:chestnut'],
                                                256)
plt.xlabel('x',fontsize=12)
plt.ylabel('y',fontsize=12)

im=plt.imshow(array2,aspect='equal',cmap = cmap)
cax2 = fig2.add_axes([1, 0.1, 0.05, 0.8])
# cax2 = fig2.add_axes([ax2.get_position().x1+0.01,ax2.get_position().y0,0.02,ax2.get_position().height])
fig2.colorbar(im, cax=cax2, shrink=2)
plt.savefig('D:\\SCIENCE\\2021\\RNF\\1\\2.png',bbox_inches = 'tight',dpi=200/fig2.get_size_inches()[1])
# -----------------------------------------------------------------------------


plt.show()
