from __future__ import division

import numpy as np
import pylab as pl
import pickle

def smooth(array):
    return [(array[i][1]+4*array[i+1][1]+6*array[i+2][1]+4*array[i+3][1]+array[i+4][1])/16 for i in range(len(array)-5)]

lambda_ = np.pi

# Calculate determinant for given values of parameters and smooth it
f = open('/home/ivan/determinant/N_1000/lambda_'+str(lambda_)+'/N_1000_lambda_'+str(lambda_)+'_T_0.45mu_0.04pi.det', 'r')
det_raw = pickle.load(f)
f.close()
det = smooth(det_raw)

mu = 0.04*np.pi
times = range(500)
T = 0.45

print det_raw

# Theoretical calculation for the comparison
teor_log_det = [np.log(T*(np.e**(1j*lambda_)-1)+1)*mu*(t+2)/(2*np.pi)-lambda_**2*np.log(np.pi*(t+2))/(4*np.pi**2) for t in times]
#teor_log_det = [np.log((T*np.e**(1j*mu*t)+1-T)/(t+1))-3 for t in times]


# Plot the real part of log(det)
pl.figure()
pl.plot(np.real(np.log(det)))
pl.plot(np.real(teor_log_det))

# Plot the imaginary part of log(det)
pl.figure()
pl.plot(np.imag(np.log(det)))
pl.plot(np.imag(teor_log_det))

# Plot the quantum corrections
pl.figure()
pl.plot(np.imag(np.log(det)-teor_log_det[:-5]))
pl.figure()
pl.plot(np.real(np.log(det)-teor_log_det[:-5]))

# Show all plots
pl.show()

