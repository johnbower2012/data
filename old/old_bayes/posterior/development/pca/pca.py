import matplotlib.pyplot as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


data1 = np.loadtxt('val.dat')
data2 = np.loadtxt('tval.dat')

#3D plot
figure = pl.figure()
ax = figure.add_subplot(111, projection='3d')

"""
x=np.random.rand(1,200)
y=np.random.rand(1,200)

x = 1+x*4
y = 1+y*4
"""

ax.scatter(data1[:,0],data1[:,1],data1[:,2],c='c')
#ax.scatter(data2[:,0],data2[:,1],data2[:,2],c='r')

pl.xlabel('x')
pl.ylabel('y')
pl.show()
