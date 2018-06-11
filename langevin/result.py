import matplotlib.pyplot as pl
import numpy as np
import math as math
from mpl_toolkits.mplot3d import Axes3D


data = np.loadtxt('result.dat')

#3D plot
figure = pl.figure()
ax = figure.add_subplot(111, projection='3d')

x=np.random.rand(1,200)
y=np.random.rand(1,200)

x = 1+x*4
y = 1+y*4

ax.scatter(data[:,0],data[:,1],data[:,2],c='k')

pl.xlabel('x')
pl.ylabel('y')


"""
#2D plot
x = np.arange(0,12,0.02)


pl.plot(x,x+np.sin(1.5*x),'k',label='f(x) = x + sin(1.5x)')
pl.plot(data2[:,0],data2[:,1],'ms',label='training points')
pl.plot(data1[:,0],data1[:,1],'r--',label='post distr mean')
pl.gca().fill_between(data1[:,0].flat,data1[:,2],data1[:,3],color="#dddddd")
pl.plot(data1[:,0],data1[:,4],'c--')
pl.plot(data1[:,0],data1[:,5],'y--')
pl.plot(data1[:,0],data1[:,6],'g--')


pl.xlabel('x')
pl.ylabel('y')
pl.grid(True)
pl.title('Posterior Functions\n\nshaded region (+/- 2 std)')
pl.legend()
#pl.xlim([0,12])
"""

pl.show()
