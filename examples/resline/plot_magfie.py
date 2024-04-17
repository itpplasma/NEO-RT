#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 17:31:12 2017

@author: ert
"""

import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('orbit_test_magfie.out')
data2 = np.loadtxt('orbit_test_magfie_neo.out')

#for k in range(1,data1.shape[1]):
for k in range(12,15):
    plt.figure()
    plt.plot(data1[:,0], data1[:,k], '-r')
    plt.plot(data2[:,0], data2[:,k], '--k')
    plt.title(k)
    
plt.show()
