#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 11:58:51 2022

@author: lukemeredith
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import minorticks_on
import functions
from functions import *


def main():
    
    m = np.array([[5,1], [3,2]])
    avg_fixation_times = []
    avg_fixation_times_A = []
    
    tau_avg = []
    tau_avg_A = []

    beta_store = np.linspace(0, 0.01 ,11)
    beta_store_2 = np.linspace(0, 0.012, 13)
   
    
    sims = 100

    pbar = tqdm(10*sims)
    
    for i, beta in enumerate(beta_store):

        fixation_time = []
        fixation_time_A = []
        
        
        for i in range(sims):
         
            t = [0] 
            n_A_store = [1]
            n_B_store = [99]
            n_A = n_A_store[0]
            n_B = n_B_store[0]
            N = n_A + n_B
    
    
            while(n_A_store[-1] != 0 and n_B_store[-1] != 0):
                
                Tplus = get_T_plus(n_A, N, m, beta)
                Tminus = get_T_minus(n_A, N, m, beta)
                
                lambda_ = Tplus + Tminus
                
                r1 = np.random.uniform(0,1)
                tau = (-1/lambda_) * np.log(r1)
                
                r2 = np.random.uniform(0,1)
                
                if r2 < Tplus/lambda_:
                    
                    n_A += 1
                    n_B -= 1
                
                else:
                    
                    n_A -= 1
                    n_B += 1
                    
                t.append(t[-1] + tau)
                n_A_store.append(n_A)
                n_B_store.append(n_B)
            
            pbar.update(1)
            
            if (n_A_store[-1] == N):
                
                fixation_time_A.append(t[-1]) 
                
                
            fixation_time.append(t[-1])
                
                
        avg_fixation_times.append(np.mean(fixation_time))
        
        if len(fixation_time_A) != 0:
            avg_fixation_times_A.append(np.mean(fixation_time_A))
        else:
            avg_fixation_times_A.append(0)
            
    pbar.close()
    
    tau_avg = np.array(avg_fixation_times)/avg_fixation_times[0]
    tau_avg_A = np.array(avg_fixation_times_A)/avg_fixation_times_A[0]

    return tau_avg, tau_avg_A


if __name__ == "__main__":
    x = main()

    tau_avg = x[0]
    tau_avg_A = x[1]
 
    np.savetxt("Tau_avg_000", tau_avg)
    np.savetxt("Tau_avg_A_000", tau_avg_A)
    print("Files saved")
    

    