# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 10:33:30 2015

@author: ivan
"""
from __future__ import division

import multiprocessing as mp
from Queue import Empty
from itertools import chain
from scipy import linalg
import numpy as np
import pylab as pl
import pickle
from time import time, sleep

PROCESSES = 4
start = time()

N = 1e5 # system size

t_RC = 3 # RC time of the detector in units of 1/e_F; 20 checked as good for 2*pi, 5 for .5*pi
dec = 8 # extend time range to allow for relaxation of detector; 8 is good

t_max = 3000  # 5000 only after clean restart!
t_min = 0
dt = 3
times = np.arange(t_min, t_max, dt)

lambda_ = -0.3*np.pi
mu = 0.15 # 2*pi/mu must be at least 12 times more than t_RC
T = 0.2


def determinant(t):
    size = t + dec*t_RC
    i,j = np.meshgrid(range(size), range(size))
    d = i - j
    tt = 1 - np.exp(-i/t_RC)  + (-1 + np.exp(-(i-t)*(np.sign(i-t) + 1)/2/t_RC))
    del i,j
    F = ((1+T*(np.exp(1j*mu*d) - 1) - np.exp(-1j*np.pi*d)+ (1-np.sign(d)**2))
         /(np.exp(2j*np.pi*d/N) - 1 + (1-np.sign(d)**2)*(1/N/(mu*T/2/np.pi + 1/2))))/N
    del d
    det = linalg.det(np.eye(size)-F+F*np.exp(1j*lambda_*tt))
    del F,tt
    return det*np.exp(-1j*lambda_*t/2) 

def work(tasks, results):
    while True:
        try:
            t = tasks.get(block=False)
            print mp.current_process().name+' gets task: ', t
            result = determinant(t)
            results.put([t, result])
            print mp.current_process().name+' puts result: ', time() - start
        except Empty: break

my_tasks = mp.Queue()
my_results = mp.Queue()

for t in times:
    my_tasks.put(t)
    
workers = [mp.Process(target=work, args=(my_tasks, my_results)) for i in range(PROCESSES)]
for each in workers:
    each.start()

rest = len(times)
dets = []
while rest:
    try:
        dets.append(my_results.get())
        rest -= 1
    except Empty :
        sleep(.01)

for each in workers:
    each.join()

dets.sort(key= lambda x : x[0])
dets = map(lambda x : x[1], dets)
f = open('/home/ivan/determinant_test/lambda_'+str(lambda_)+'T_'+str(T)+'mu_'+str(mu)+'.det', 'w')
pickle.dump(dets, f)
f.close()

dets = np.log(dets)

pl.plot(np.real(dets))
pl.plot(-2*(lambda_/(2*np.pi))**2*np.log(1+times/t_RC))
pl.figure()
pl.plot(np.imag(dets))

pl.show()