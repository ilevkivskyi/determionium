from __future__ import division

import multiprocessing as mp
from Queue import Empty
from time import time, sleep
from scipy import linalg
import numpy as np
import pylab as pl
import pickle
import os

# Defining auxiliary functions
def delta(x1, x2):
    if x1 == x2: return 1
    else : return 0

def theta(x):
    if x > 0 : return 1
    else : return 0

def ft(n, t):
    if n != 0 :
        return (1/N) * (np.e**(1j*2*np.pi*n*t/N)-1)/(np.e**(1j*n*2*np.pi/N)-1)
    else :
        return t / N

# Main function
def determ(t, T, lambda_):
    
    S = [[(1-T)        , np.sqrt(T*(1-T))],
         [np.sqrt(T*(1-T)),            T]]
    
    L_t = [[ft(m-n, t)*theta(N/2 - abs(m-n)) for n in range(N)] for m in range(N)]
    L_texp = linalg.expm(1j*lambda_*np.kron(L_t,S))
    mat = np.eye(2*N) - F + np.matrix(L_texp)*F
     
    return np.linalg.det(np.array(mat))*np.e**(-1j*lambda_*t/2)

def work(tasks, results):
    while True:
        try:
            args = tasks.get(block=False)
            print mp.current_process().name + ' get task: ', time() - timer0
            result = determ(*args)
            results.put([args[0], result])
            print mp.current_process().name + ' put result: ', time() - timer0
        except Empty: break

if __name__ == '__main__':
    PROCESSES = 5
    
    # Parameters of the problem
    N = 1000
    mu = 0.04*np.pi
    t_values = range(500)
    T_values = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.475, 0.49, 0.495, 0.5, 0.505, 0.51, 0.525, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]
    lambda_values = [np.pi/3, np.pi, 2*np.pi]
    
    # Constructing distribution function matrix F once for all
    f_0 = [[delta(m,n)*theta(N/2-n) for n in range(N)] for m in range(N)]
    f_mu = [[delta(m,n)*theta(N/2+mu*N/(2*np.pi)-n) for n in range(N)] for m in range(N)]
    L = [[ 0, 0],
         [ 0, 1]]
    U = [[ 1, 0],
         [ 0, 0]]
    F = np.matrix(np.kron(f_0, U)+np.kron(f_mu, L))
    
    my_tasks = mp.Queue()
    my_results = mp.Queue()
    
    try:
        os.mkdir('/home/ivan/determinant/N_'+str(N))
    except OSError : pass
    for lambda_ in lambda_values:
        try:
            os.mkdir('/home/ivan/determinant/N_'+str(N)+'/lambda_'+str(lambda_))
        except OSError : pass
        for T in T_values:
            # Doing calculations for all times...
            for t in t_values:
                my_tasks.put((t, T, lambda_))
            workers = [mp.Process(target=work, args=(my_tasks, my_results)) for i in range(PROCESSES)]
            timer0 = time()
            for each in workers:
                each.start()
            rest = len(t_values)
            det_raw = []
            while rest:
                try:
                    det_raw.append(my_results.get())
                    rest -= 1
                except Empty : pass
            for each in workers:
                each.join()
            # ... sort and save the results in a file
            det_raw.sort(key= lambda x : x[0])
            f = open('/home/ivan/determinant/N_'+str(N)+'/lambda_'+str(lambda_)+'/N_'+str(N)+'_lambda_'+str(lambda_)+'_T_'+str(T)+'mu_0.04pi'+'.det', 'w')
            pickle.dump(det_raw, f)
            f.close()
