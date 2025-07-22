import numpy as np
import math, cmath

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as plticker
from matplotlib import cm

from scipy.spatial import Delaunay
import surf2stl_python.surf2stl as s2s

N = 9
a = 0.4
row, col = 30, 30
writeSTL = False

def calcZ1(x, y, k, n):
    return cmath.exp(1j*(2*cmath.pi*k/n)) * (cmath.cos(x+y*1j)**(2/n))

def calcZ2(x, y, k, n):
    return cmath.exp(1j*(2*cmath.pi*k/n)) * (cmath.sin(x+y*1j)**(2/n))

def calcZ1Real(x, y, k, n):
    return (calcZ1(x, y, k, n)).real

def calcZ2Real(x, y, k, n):
    return (calcZ2(x, y, k, n)).real

def calcZ(x, y, k1_, k2_, n, a_):
    z1 = calcZ1(x, y, k1, n)
    z2 = calcZ2(x, y, k2, n)
    return z1.imag * math.cos(a_) + z2.imag*math.sin(a_)

# set param range
x = np.linspace(0, math.pi/2, col)
y = np.linspace(-math.pi/2, math.pi/2, row)
x, y = np.meshgrid(x, y)

# init graph
fig = plt.figure(figsize=(18,8))

for n in range(2, N):
    ax = fig.add_subplot(2, 4, n - 1, projection='3d')
    ax.view_init(elev=15, azim=15)
    ax.set_title("n=%d" % n)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_locator(loc)
    ax.zaxis.set_major_locator(loc)

    count = 0
    for k1 in range(n):
        for k2 in range(n):
            # calc X, Y, Z values
            X = np.frompyfunc(calcZ1Real, 4, 1)(x, y, k1, n).astype('float32')
            Y = np.frompyfunc(calcZ2Real, 4, 1)(x, y, k2, n).astype('float32')
            Z = np.frompyfunc(calcZ, 6, 1)(x, y, k1, k2, n, a).astype('float32')

            ax.plot_surface(X, Y, Z, cmap=cm.ocean, linewidth=0)

            # write to a STL file
            if writeSTL:
                X_ = X.flatten()
                Y_ = Y.flatten()
                Z_ = Z.flatten()
                delaunay_tri = Delaunay(np.array([x_, y_]).T)
                s2s.tri_write('output_calabi-yau_n%d_%d.stl'% (n, count), X_, Y_, Z_, delaunay_tri)

            count += 1
