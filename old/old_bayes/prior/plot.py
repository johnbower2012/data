import matplotlib.pyplot as pl
import numpy as np

data = np.loadtxt('prior.dat')

pl.plot(data[:,0],data[:,1], 'ro')
pl.plot(data[:,0],data[:,2], 'mo')
pl.plot(data[:,0],data[:,3], 'go')
pl.plot(data[:,0],data[:,4], 'yo')
pl.plot(data[:,0],data[:,5], 'bo')

pl.xlabel('x')
pl.ylabel('f(x)')
pl.grid(True)
pl.title('Prior Functions \nmean=0, s=1')
pl.legend()
pl.show()


