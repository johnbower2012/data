import matplotlib.pyplot as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


data1 = np.loadtxt('hyperp.dat')
"""
#3D plot
figure = pl.figure()
ax = figure.add_subplot(111, projection='3d')

x=np.random.rand(1,200)
y=np.random.rand(1,200)

x = 1+x*4
y = 1+y*4

ax.scatter(data2[:,0],data2[:,1],data2[:,2],c='r')
ax.scatter(data1[:,0],data1[:,1],data1[:,2],c='c')
ax.scatter(data1[:,0],data1[:,1],data1[:,5],c='g')
#ax.scatter(data1[:,0],data1[:,1],data1[:,6],c='y')
#ax.scatter(data1[:,0],data1[:,1],data1[:,7],c='m')

pl.xlabel('x')
pl.ylabel('y')

"""
#2D plot
x = np.arange(0,5,0.02)


pl.plot(data1[:,0],-data1[:,1],'co')
pl.plot(data1[:,0],data1[:,2],'ro')



pl.xlabel('x')
pl.ylabel('y')
#pl.yscale('log')
pl.grid(True)
pl.legend()
#pl.xlim([0,12])

pl.show()
