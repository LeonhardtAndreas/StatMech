#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import gamma


x = np.arange(21)
y = np.zeros_like(x,dtype=float)
z = np.zeros_like(x,dtype=float)
R = 1 
for i,n in enumerate(x):
    y[i] = np.sqrt(np.pi)**n / gamma(n/2+1) 
    z[i] = np.sqrt(np.pi)**n / gamma(n/2+1)/2**n 


fig, [ax1,ax2] = plt.subplots(2,1)
ax1.plot(x,y,'o')
ax1.set_xlabel(r'$n$')
ax1.set_ylabel(r'$V_n(1)$')
ax1.set_xlim([0,20])
ax1.set_xticks(x)

ax2.plot(x,z,'o')
ax2.set_xlabel(r'$n$')
ax2.set_ylabel(r'$V_n(\frac{1}{2})$')
ax2.set_xlim([0,20])
ax2.set_xticks(x)
plt.show()
