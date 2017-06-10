# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 10:59:39 2015

@author: ivan
"""
from __future__ import division

import numpy as np
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

title = 'rc_det_test'

f_out = title + '.dets'
f_aux = title + '.aux'

with open(f_aux) as f:
    t_vals = map(float, f.readline().split())
    t_RC_vals = map(float, f.readline().split())
    la_vals = map(float, f.readline().split())
    mu_vals = map(float, f.readline().split())
    T_vals = map(float, f.readline().split())

with open(f_out) as f:
    lines = [line.split() for line in f.readlines()]

lines = [map(float, line) for line in lines]
lines.sort(key = lambda x:x[0])
data0 = np.array(lines)

data = data0[np.where(np.abs(data0.T[2]-2*np.pi)<0.001)] # make crossections

dets = []
for t in t_vals:
    row = list(data[np.where(data.T[0]==t)])
    row.sort(key = lambda x:x[1]) # sort by second param. for each time
    row = map(lambda x:x[-2]+1j*x[-1], row)  # -2 for real, -1 for imaginary
    dets.append(row)
dets = np.array(dets)


X, Y = np.meshgrid(t_vals, t_RC_vals)

dets = np.log(np.exp(dets))

redets = np.real(dets)
imdets = np.imag(dets)

fig = pl.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, redets.T, rstride=1, cstride=2, 
                       cmap=cm.coolwarm, linewidth=0, antialiased=False)

ax.set_zlim(redets.min()-0.1, redets.max()+0.1)
fig.colorbar(surf, shrink=0.5, aspect=8)

pl.figure()
pl.imshow(redets.T[22:,100:], aspect = 20)

pl.show()