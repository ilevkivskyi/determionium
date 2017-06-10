# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:51:56 2015

@author: ivan
"""
from __future__ import division

from itertools import product
import numpy as np

title = '/home/levkivsk/c_det_test'

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