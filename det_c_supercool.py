# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 21:26:28 2015

@author: ivan
"""
from __future__ import division

import numpy as np
import pylab as pl
import os
from itertools import product
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

title = '/home/ivan/c_det_test'

f_in = title + '.tasks'
f_out = title + '.dets'
f_aux = title + '.aux'

t_vals = np.arange(300)
t_RC_vals = [3,8,20]
la_vals = np.arange(0.1,2, 0.1)*np.pi
mu_vals = [0.15]
T_vals = [0.0]

with open(f_aux, 'w') as f:
    for meta in [t_vals, t_RC_vals, la_vals, mu_vals, T_vals]:
        f.write(' '.join(map(str, meta))+'\n')

with open(f_in, 'w') as f:
    for task in product(t_vals, t_RC_vals, la_vals, mu_vals, T_vals):
        f.write(' '.join(map(str, task))+'\n')

os.system('mpiexec -n 5 /home/ivan/det_return '+f_in+' '+f_out)

with open(f_out) as f:
    lines = [line.split() for line in f.readlines()]

lines = [map(float, line) for line in lines]
lines.sort(key = lambda x:x[0])
data = np.array(lines)

data = data[np.where(data.T[1]==3)] # make a crossection

dets = []
for t in t_vals:
    row = list(data[np.where(data.T[0]==t)])
    row.sort(key = lambda x:x[2]) # sort by second param. for each time
    row = map(lambda x:x[-2], row)  # -2 for real, -1 for imaginary
    dets.append(row)
dets = np.array(dets)

X, Y = np.meshgrid(t_vals, la_vals)

fig = pl.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, dets.T, rstride=1, cstride=2, 
                       cmap=cm.coolwarm, linewidth=0, antialiased=False)

ax.set_zlim(dets.min()-0.1, dets.max()+0.1)
fig.colorbar(surf, shrink=0.5, aspect=8)

pl.figure()
pl.imshow(dets)

pl.show()