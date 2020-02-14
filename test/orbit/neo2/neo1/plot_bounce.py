#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 17:31:12 2017

@author: ert
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('standalone2/test_bounce.dat')

t = np.linspace(0.0,1.0,len(data))

plt.plot(np.sqrt(data[:,0])*np.cos(data[:,2]), np.sqrt(data[:,0])*np.sin(data[:,2]))
plt.plot(np.sqrt(data[:,4])*np.cos(data[:,6]), np.sqrt(data[:,4])*np.sin(data[:,6]))
plt.plot(np.sqrt(data[:,8])*np.cos(data[:,10]), np.sqrt(data[:,8])*np.sin(data[:,10]), '--')
plt.xlim(-1,1)
plt.ylim(-1,1)

