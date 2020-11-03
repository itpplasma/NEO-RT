from numpy import *
from matplotlib.pyplot import *

for k in range(7, 11+1):
  data = loadtxt(f'fort.{1000+k}')
  plot(data[:, 0], data[:, 2], '.')
