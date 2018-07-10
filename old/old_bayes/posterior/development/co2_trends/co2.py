import matplotlib.pyplot as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


test_data = np.loadtxt('co2_out.dat')
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

pl.plot(train_data[:,0],train_data[:,1],'bo',label='training points')

pl.plot(test_data[:,0],test_data[:,1],'rs',label='mean of distribution')
pl.gca().fill_between(test_data[:,0].flat,test_data[:,2],test_data[:,3],color="#dddddd", label='2 sigma')
pl.plot(test_data[:,0],test_data[:,4],'co') 
pl.plot(test_data[:,0],test_data[:,5],'mo') 
pl.plot(test_data[:,0],test_data[:,6],'yo')

pl.xlim([1955,2030])
pl.ylim([300,460])
pl.xlabel('Year')
pl.ylabel('CO2 ppm')
pl.grid(True)
pl.legend(loc='upper left')
pl.title('Atmospheric CO2 in ppm')

pl.show()
