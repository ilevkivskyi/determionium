from __future__ import division
import numpy as np
from scipy.optimize import fsolve

EPSILON = 0.001
SEARCH_ZERO = 2
SEARCH_SECOND_MAX = 10

#T_list = [0.15,0.2,0.25,0.3,0.35,0.4,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]
T_list = [0.15,0.2,0.25,0.3,0.35,0.4]
T_list = [0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]

# derivative operator
def D(f):
    def derivative(x, *args):
        return (f(x + EPSILON, *args) - f(x, *args)) / EPSILON
    return derivative

# Visibility asymptotics for Gaussian and non-Gaussian cases
# bias mu is measured in units of 2v/L
def V(mu, T, stat):
    if stat == 'G':
        return (T*np.cos(T*mu) + (T - 1)*T*np.pi*np.sin(T*mu)/2) * np.e**((T - 1)*T*np.pi*mu/2)

# non-Gaussian case two asymptotics, above and belowe phase transition
    if stat == 'nG':
        if T > 0.5: 
            return (np.log(2*T - 1)*np.sin(mu)/np.pi + np.cos(mu)) * np.e**(mu*np.log(2*T - 1)/np.pi)
        else: 
            return (1 + mu*np.log(1 - 2*T)/np.pi) * np.e**(mu*np.log(1 - 2*T)/np.pi)

first_zero = {'G' : [], 'nG' : []}
second_max_pos = {'G' : [], 'nG' : []}
second_max_value = {'G' : [], 'nG' : []}

for stat in ['G', 'nG']:
    for i,T in enumerate(T_list):
        first_zero[stat].append(fsolve(V, SEARCH_ZERO, (T, stat))[0])
        second_max_pos[stat].append(fsolve(D(V), SEARCH_SECOND_MAX, (T, stat))[0])
        second_max_value[stat].append(V(second_max_pos[stat][i], T, stat))

print first_zero['G']
print first_zero['nG']

print second_max_value['G']
print second_max_value['nG']
