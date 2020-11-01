from numpy import *
from matplotlib.pyplot import *

data = loadtxt('fort.999')
plot(data[:, 0], data[:, 2])

print(mean(data[:,1]))
