import matplotlib.pyplot as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt('delY.dat')
data1 = np.loadtxt('I211_J321.dat')

pl.plot(data,data1[:,0],'c-',label='l0=24.84')
pl.plot(data,data1[:,1],'r-',label='l1=3.23')
pl.plot(data,data1[:,2],'m-',label='l2=1.01')
pl.plot(data,data1[:,2],'y-',label='l3=0.92')



pl.legend()
pl.grid()
pl.title('I211_J321.dat')
pl.xlabel('delY')
pl.ylabel('coeff')
pl.show()
