# -*- coding: utf-8 -*-
"""
Created on Sun Feb 22 14:59:08 2015

@author: ivan
"""
from __future__ import division

import numpy as np
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

title5 = '/home/ivan/e_det_test'
title = '/home/ivan/d_det_test'

f_out = title + '.dets'
f_aux = title + '.aux'
f_out5 = title5 + '.dets'

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
data = np.array(lines)

data0 = data[np.where(data.T[4]==0.0)] # make crossections


with open(f_out5) as f:
    lines = [line.split() for line in f.readlines()]

lines = [map(float, line) for line in lines]
lines.sort(key = lambda x:x[0])
data5 = np.array(lines)


data2 = data5[np.where(data5.T[4]==0.5)] 

dets0 = []
for t in t_vals:
    row = list(data0[np.where(data0.T[0]==t)])
    row.sort(key = lambda x:x[2]) # sort by second param. for each time
    row = map(lambda x:x[-2]+1j*x[-1], row)  # -2 for real, -1 for imaginary
    dets0.append(row)
dets0 = np.array(dets0)

dets2 = []
for t in t_vals:
    row = list(data2[np.where(data2.T[0]==t)])
    row.sort(key = lambda x:x[2])
    row = map(lambda x:x[-2]+1j*x[-1], row)
    dets2.append(row)
dets2 = np.array(dets2)

dets  = dets2 - dets0
edets = np.exp(dets)
f = np.fft
fts = f.fftshift(f.hfft(edets, axis=0), axes=0)

for i in range(len(la_vals)):
    fts[:,i] = fts[:,i]/np.cumsum(fts[:,i])[-1]
    fts[:,i] = fts[:,i]/max(fts[:,i])

e_vals = np.linspace(-np.pi,np.pi, 2*len(t_vals)-2)

X, Y = np.meshgrid(e_vals, la_vals)

redets = np.real(fts)

fig = pl.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, redets.T, rstride=1, cstride=2, 
                       cmap=cm.coolwarm, linewidth=0, antialiased=False)

ax.set_zlim(redets.min()-0.1, redets.max()+0.1)
fig.colorbar(surf, shrink=0.5, aspect=8)

pl.figure()
pl.imshow(redets.T, aspect = 5)

pl.show()