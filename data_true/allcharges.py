import matplotlib.pyplot as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt('delY.dat')
data1 = np.loadtxt('allcharges.dat')

pl.plot(data,data1[:,0],'c-',label='l0=10.44')
pl.plot(data,data1[:,1],'r-',label='l1=3.40')
pl.plot(data,data1[:,2],'m-',label='l2=0.34')
pl.plot(data,data1[:,2],'y-',label='l3=0.29')

pl.legend()
pl.grid()
pl.title('allcharges.dat')
pl.xlabel('delY')
pl.ylabel('coeff')
pl.show()
