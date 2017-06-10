from __future__ import division
from itertools import chain
from pylab import plot, show, figure, contourf
import numpy as np
import pickle

NN = 2*(2*(2*(3072-1)-1)-1)

def smooth(array):
    return np.array([(array[i]+4*array[i+1]+6*array[i+2]+4*array[i+3]+array[i+4])/16 for i in range(len(array)-5)])

def interpolate_2(arr):
    newarr = []
    for i in range(len(arr) - 1):
        newarr.append(arr[i])
        newarr.append((arr[i] + arr[i+1]) / 2)
    return np.array(newarr)

def interp(n, arr):
    for i in range(n):
        arr = interpolate_2(arr)
    return arr

def gaus(x):
    s = np.linspace(0,1,1000)[1:]
    return sum((1-np.cos(x*s))*(1-s)/s**2 for s in s)/1000

eta = .27
mu = 0.015*np.pi
time = np.arange(0, 3072)

plot(gaus(time))
figure()

T_values = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3,
                    0.35, 0.4, 0.45, 0.475, 0.5, 0.7, 1.0])


home = '/home/ivan/Pandemonium/out_4000_pi_0015_d_'
with open(home + '1.txt') as f1, \
     open(home + '2.txt') as f2, \
     open(home + '3.txt') as f3:
    log_correlators = \
    [[complex(num.strip(',').replace('+-', '-')) for num in line.split()]
                                                 for line in chain(f1, f2, f3)]

correlators = np.exp(np.array(log_correlators).T[:-1])
#correlators = [interp(3, c) for c in correlators]
correlators_n = [list(reversed(np.conj(corr / correlators[0]))) + list(corr / correlators[0]) for corr in correlators]

T_ext = np.array([0.025, 0.075, 0.525, 0.55, 0.6, 0.65, 0.7,
                    0.75, 0.8, 0.85, 0.9, 0.95])

home = '/home/ivan/Pandemonium/out_4000_pi_0015_d_ext_'
with open(home + '1.txt') as f1, \
     open(home + '2.txt') as f2, \
     open(home + '3.txt') as f3:
    log_correlators = \
    [[complex(num.strip(',').replace('+-', '-')) for num in line.split()]
                                                 for line in chain(f1, f2, f3)]

correlators_ext = np.exp(np.array(log_correlators).T[:-1])
correlators_n_ext = [list(reversed(np.conj(corr / correlators_ext[0]))) + list(corr / correlators_ext[0]) for corr in correlators_ext]


c = correlators_n[11]

T = 0.9
time2 = np.arange(-3072, 3072)
#c = np.exp(np.log(1-2*T+0j)*mu*abs(time2)/(2*np.pi))
NN = 3072

c = np.exp(1j*T*mu*time2-T*(1-T)*gaus(mu*time2))

plot(np.abs(smooth(np.diff(np.abs(np.imag(c)))))[3070:3470])
figure()

vis = [sum(m*c[t]*c[t-m]*(  1/((t -NN - eta*1j)*(t-NN-m - eta*1j)) - 1/((t -NN + eta*1j)*(t-NN-m + eta*1j)) ) for t in chain(range(m, NN+1), range(NN+1, NN+m+1), range(NN+m+1, 2*NN) ) ) for m in range(1, 500)]
#vis2 = correlators/correlators[0] - np.conj(correlators/correlators[0])

print vis

plot(np.abs(np.diff(vis)-0*vis[0]))

#save_data = np.abs(smooth(np.diff(np.abs(np.imag(c)))))[3070:3470]
save_data = np.abs(np.diff(vis)-0*vis[0])

f = open('/home/ivan/gaus_04.txt','w')
for d in save_data:
    f.write(str(d)+', ')
f.close()

Ts = [0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3,
                    0.35, 0.4, 0.45, 0.475, 0.5, 0.525, 0.55, 0.6, 0.65, 0.7,
                    0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
cs = correlators_n[0:1] + correlators_n_ext[0:1] + correlators_n[1:2] + correlators_n_ext[1:2] + correlators_n[2:12] + correlators_n_ext[2:] + correlators_n[13:]

#sec_zero = [nan, nan, nan, nan, nan, nan, nan, nan, nan, 248.5, 236, 225, 219, 214.5, 211, 208.5, 206.5, 204.5, 202.5, 199]

#thi_zero = [nan, nan, nan, nan, nan, nan, nan, nan, nan, 391, 378, 364.5, 356.5, 351.5, 347, 344, 341, 339, 336.5, 332.5]

#res = []
#for i in range(len(cs)):
    #c = cs[i]
    #if i != 13:
        #vis = np.abs(np.diff([sum(m*c[t]*c[t-m]*(  1/((t -NN - eta*1j)*(t-NN-m - eta*1j)) - 1/((t -NN + eta*1j)*(t-NN-m + eta*1j)) ) 
        #for t in chain(range(m, NN+1), range(NN+1, NN+m+1), range(NN+m+1, 2*NN) ) ) for m in range(1, 500)]))
    #else:
        #vis = np.abs(smooth(np.diff(np.abs(np.imag(c)))))[3070:3570]
    #res.append(vis/vis[4])

#f = open('/home/ivan/cool','w')
#pickle.dump(res,f)
#f.close()

f = open('/home/ivan/cool','r')
cool = pickle.load(f)
f.close()

cool[13] = cool[13][:-2]

cool = cool[1:][:]

for c in cool:
    print len(c), '\n'

figure()
contourf(cool, np.linspace(0,1.1,20))

#figure()
#plot(cool[6])

#plot(np.abs(8*np.array(vis2[4])))

show()
