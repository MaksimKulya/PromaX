from matplotlib import pyplot as plt
import numpy as np

def peak(x, c):
    return np.exp(-np.power(x - c, 2) / 16.0)

def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def half_max_x(x, y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[0], half),
            lin_interp(x, y, zero_crossings_i[1], half)]

# make some fake data
x=np.linspace(0,20,50)
y=peak(x,10)

# find the two crossing points
hmx = half_max_x(x,y)

# print the answer
fwhm = hmx[1] - hmx[0]
print("FWHM:{:.3f}".format(fwhm))

# a convincing plot
half = max(y)/2.0
plt.plot(x,y)
plt.plot(hmx, [half, half])

plt.show()

