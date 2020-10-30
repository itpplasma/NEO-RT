from numpy import *
from matplotlib.pyplot import *

data = loadtxt('fort.99')
plot(data[:, 1], data[:, 3])
