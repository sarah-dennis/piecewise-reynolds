# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:36:17 2024

@author: sarah
"""
import numpy as np
import readwrite as rw
import stokes_pressure as pressure
#-------------------------------------------------------------------------------
def stokes_cnvg_self(Ex, args, U, Q, Re, N_min, Ns, N_max, p_err=True):
    ex_min = Ex(args, U, Q, Re, N_min)
    ex_max = Ex(args, U, Q, Re, N_max)
    size=ex_min.Ny*ex_min.Nx
    # Load max grid for 'true' values
    u_max, v_max, psi_max, past_iters = rw.read_stokes(ex_max.filestr+".csv", ex_max.Nx * ex_max.Ny)
    psi_max = psi_max.reshape((ex_max.Ny,ex_max.Nx))


    if p_err:
        p_max = pressure.pressure(ex_max, u_max, v_max)
        dp_max = pressure.get_dp(ex_max, p_max) 
        p_max = p_max.reshape((ex_max.Ny,ex_max.Nx))


        
    mult_max = int(N_max/N_min)
    
    err_psi = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
    err_psi_inf = np.zeros(len(Ns)+1)
    err_psi_l1 = np.zeros(len(Ns)+1)
    err_psi_l2 = np.zeros(len(Ns)+1)
    

    if p_err:
        err_p = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
        err_p_inf = np.zeros(len(Ns)+1)
        err_p_l1 = np.zeros(len(Ns)+1)
        err_p_l2 = np.zeros(len(Ns)+1)
        err_dp = np.zeros(len(Ns)+1)
    
    for n in range(len(Ns)+1):
        if n == 0:
            ex_n = ex_min
            mult = 1
        else:
            ex_n = Ex(args, U, Q, Re, Ns[n-1])
            mult = int(Ns[n-1]/N_min)
            
            
        u_n, v_n, psi_n, past_iters = rw.read_stokes(ex_n.filestr+".csv", ex_n.Nx * ex_n.Ny)
        psi_n=psi_n.reshape((ex_n.Ny, ex_n.Nx))

        if p_err:
            p_n = pressure.pressure(ex_n, u_n, v_n)
            dp_n = pressure.get_dp(ex_n, p_n)
            p_n=p_n.reshape((ex_n.Ny,ex_n.Nx))
            
        for k_min in range(size):
        # all indices (i,j) on grid N_min
            i = k_min % ex_min.Nx
            j = (k_min // ex_min.Nx)
            
            
            
            # -> indices on grid N || N_min
            i_n = mult * i
            j_n = mult * j

            # -> indices on grid N_max || N_min
            i_max = mult_max * i
            j_max = mult_max * j
            if ex_max.space[j_max,i_max] != 1 or ex_n.space[j, i] !=1:
                continue
            err_psi_n_ij = abs(psi_max[j_max,i_max] - psi_n[j_n,i_n])
            err_psi[n,j,i] = err_psi_n_ij
            err_psi_l1[n] += err_psi_n_ij
            err_psi_l2[n] += err_psi_n_ij**2
            
            if p_err:
                err_p[n,j,i] = abs(p_max[j_max,i_max] - p_n[j_n,i_n])
                err_p_l1[n]+= err_p[n,j,i] 
                err_p_l2[n]+= (err_p[n,j,i]**2) 
        
        err_psi_inf[n] = np.max(err_psi[n])
        err_psi_l1[n] /= size
        err_psi_l2[n] = np.sqrt(err_psi_l2[n]/size)
        
        if p_err:
            err_p_inf[n] = np.max(err_p[n])
            err_p_l1[n] /= size
            err_p_l2[n] = np.sqrt(err_p_l2[n]/size)
            err_dp[n]= abs(dp_max-dp_n)
            
    l1_rate = convg_rate(err_psi_l1)
    l2_rate = convg_rate(err_psi_l2)
    inf_rate = convg_rate(err_psi_inf)
    cnvg_rates = np.stack([l1_rate, l2_rate, inf_rate], axis=0)
    
    if p_err:
        p_l1_rate = convg_rate(err_p_l1)
        p_l2_rate = convg_rate(err_p_l2)
        p_inf_rate = convg_rate(err_p_inf)
        dp_rate = convg_rate(err_dp)
        p_cnvg_rates = np.stack([p_l1_rate, p_l2_rate, p_inf_rate, dp_rate], axis=0)
        
    print("stream cnvg rates")
    print("l1: " + np.array2string(cnvg_rates[0], precision=2))
    print("l2: " + np.array2string(cnvg_rates[1], precision=2))
    print("linf" + np.array2string(cnvg_rates[2], precision=2))
    
    if p_err:
        print("pressure cnvg rates")
        print("l1: " + np.array2string(p_cnvg_rates[0], precision=2))
        print("l2: " + np.array2string(p_cnvg_rates[1], precision=2))
        print("linf" + np.array2string(p_cnvg_rates[2], precision=2))
        print("dp: " + np.array2string(p_cnvg_rates[3], precision=2))

    return err_psi_l1, err_psi_l2, err_psi_inf, cnvg_rates, ex_min

#------------------------------------------------------------------------------
    
def convg_rate(errs):
    n = len(errs)
    rates = np.zeros(n-1)
    
    for k in range(n-1):
        rates[k]=errs[k+1]/errs[k]
    
    return rates

#------------------------------------------------------------------------------


