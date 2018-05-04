import matplotlib.pyplot as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


test_data = np.loadtxt('testset.dat')
train_data = np.loadtxt('trainset.dat')

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


pl.plot(train_data[:,0],train_data[:,4],'bo',label='training points')
pl.plot(test_data[:,0],test_data[:,4],'r--',label='mean of distribution')
pl.gca().fill_between(test_data[:,0].flat,test_data[:,5],test_data[:,6],color="#dddddd", label='2 sigma')
pl.plot(test_data[:,0],test_data[:,7],'c-') 


pl.xlabel('x')
pl.ylabel('y')
pl.grid(True)
pl.legend(loc='upper right')
pl.title('star_KK Balance Function')

pl.show()
