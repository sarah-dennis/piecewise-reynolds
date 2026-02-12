
# -*- coding: utf-8 -*-
"""
Created on Mon May 19 17:32:18 2025

@author: sarah
"""
import numpy as np
import time
from multiprocessing import shared_memory, Pool
import reyn_boundary as bc
def make_rhs(height, BC): 
    N = height.N_regions 
    rhs = np.zeros(2*N-1)
    
    if isinstance(BC, bc.Fixed):
        rhs[0] = -BC.p0/height.widths[0] # = dp{0} - p{1}/dx{0}
        
    elif isinstance(BC, bc.Mixed):
        h0 = height.h_steps[0]
        rhs[0] = -12*BC.Q*h0**-3 + 6*BC.U*h0**-2 # = dp{0}
        

    rhs[N-1] = BC.pN/height.widths[-1] # = dp{N-1} + p{N-1}/dx{N-1}
    
    sixU = 6*BC.U
    h_steps = height.h_steps
   
    for k in range (N-1):
        rhs[N + k] = (h_steps[k+1] - h_steps[k]) * sixU
       
    return rhs

def schur_solve(height, BC):

    N = height.N_regions

    rhs = make_rhs(height, BC)

    S, S_prod = get_S(height, BC)
    D = get_D(height, BC, S)
    
    h_steps = height.h_steps
    widths = height.widths
    
    p_peaks = np.zeros(N-1) # interior peaks only
    
    for i in range(N-1):
        p_peak_ij = 0

        j=0
        k_inv_ij = K_inv_ij(D, S_prod, i, j)
        p_peak_i_j = k_inv_ij * (rhs[N+j] + rhs[j] *  h_steps[j]**3)
        p_peak_ij += p_peak_i_j
        
        j=N-2
        k_inv_ij = K_inv_ij(D, S_prod, i, j)
        p_peak_i_j = k_inv_ij * (rhs[N+j] - rhs[j+1] * h_steps[j+1]**3)
        p_peak_ij += p_peak_i_j
        
        for j in range(1,N-2):
            
            if rhs[N+j]!=0:
                k_inv_ij = K_inv_ij(D, S_prod, i, j)
     
                p_peak_i_j = k_inv_ij * rhs[N+j]
  
                p_peak_ij += p_peak_i_j
                
        p_peaks[i] = p_peak_ij
        

    p_slopes = np.zeros(N)

    if isinstance(BC, bc.Fixed):
        p0 = BC.p0
        p_slopes[0] = (p_peaks[0] - p0)/widths[0]
        
    elif isinstance(BC, bc.Mixed):
        p0 = p_peaks[0]  - rhs[0] * widths[0]
        p_slopes[0] = rhs[0]
    
    for i in range(1, N-1):
        p_slopes[i] = (p_peaks[i] - p_peaks[i-1])/widths[i]
        
    p_slopes[N-1] = (BC.pN - p_peaks[N-2])/widths[N-1]

    # print('schur time: ', tf-t0)
    ps = make_ps(height, BC, p_slopes, p_peaks)
    # tF = time.time()
    # print('schur total time:', tF-t0)
    return ps# , tf-t0


def schur_solve_parallel(height, BC):

    N = height.N_regions

    rhs = make_rhs(height, BC)
   
    S, S_prod = get_S(height, BC)
    D = get_D(height, BC, S)
    
    h_steps = height.h_steps
    widths = height.widths
    
    p_peaks = np.zeros(N-1) # interior peaks only
    
    smem_h_steps, smem_rhs, smem_D, smem_S_prod = create_shared_memory(h_steps, rhs, D, S_prod)
    
    with Pool(processes=min(12,int(height.N**(1/2)))) as pool:
        p_peaks = pool.starmap(outer_process, ((i, N, smem_h_steps, smem_rhs, smem_D, smem_S_prod) for i in range(N-1))) 
        
        # for i in range(N-1):    
        #     p_peak_i_js = pool.starmap(inner_process, ((i, j, N, smem_h_steps, smem_rhs, smem_D, smem_S_prod) for j in range(N-1)))
        #     p_peaks[i] = sum(p_peak_i_js)
            
            
    smem_h_steps.shm.close()
    smem_h_steps.shm.unlink()
    smem_rhs.shm.close() 
    smem_rhs.shm.unlink()
    smem_D.shm.close()
    smem_D.shm.unlink()
    smem_S_prod.shm.close()
    smem_S_prod.shm.unlink()
    
    p_slopes = np.zeros(N)

    if isinstance(BC, bc.Fixed):
        p0 = BC.p0
        p_slopes[0] = (p_peaks[0] - p0)/widths[0]
        
    elif isinstance(BC, bc.Mixed):
        p0 = p_peaks[0]  - rhs[0] * widths[0]
        p_slopes[0] = rhs[0]
    
    for i in range(1, N-1):
        p_slopes[i] = (p_peaks[i] - p_peaks[i-1])/widths[i]
        
    p_slopes[N-1] = (BC.pN - p_peaks[N-2])/widths[N-1]

    # print('schur time: ', tf-t0)
    ps = make_ps(height, BC, p_slopes, p_peaks)
    # tF = time.time()
    # print('schur total time:', tF-t0)
    return ps# , tf-t0

def create_shared_memory(h_steps, rhs, D, S_prod):
    smem_h_steps = shared_memory.ShareableList(h_steps.tolist(),name='smem_h_steps')
    smem_rhs = shared_memory.ShareableList(rhs.tolist(),name='smem_rhs')
    smem_D = shared_memory.ShareableList(D.tolist(),name='smem_D')
    smem_S_prod = shared_memory.ShareableList(S_prod.tolist(),name='smem_S_prod')

    return smem_h_steps, smem_rhs, smem_D, smem_S_prod

def inner_process(i, j, N, h_steps, rhs, D, S_prod):
    k_inv_ij = K_inv_ij(D, S_prod, i, j)
    
    if j == 0:
        p_peak_ij = k_inv_ij * (rhs[N] + rhs[0] *  h_steps[0]**3)
      
    elif j == N-2:
        p_peak_ij = k_inv_ij * (rhs[2*N-2] - rhs[N-1] * h_steps[N-1]**3)
        
    else: 
        p_peak_ij = k_inv_ij * rhs[N+j]
    return p_peak_ij

def outer_process(i, N, h_steps, rhs, D, S_prod):
    p_peak_i_j =0
    for j in range(N-1):
        p_peak_i_j += inner_process(i, j, N, h_steps, rhs, D, S_prod)
    return p_peak_i_j



def K_inv_ij(D, S_prod, i, j):
    if i == j:
        k_inv_ij = D[i]
        
    elif i < j:
        k_inv_ij = D[i] * (S_prod[j-1]/S_prod[i-1]) 
        
    else:
        k_inv_ij = D[j] * (S_prod[i-1]/S_prod[j-1])
    return k_inv_ij

# schur complement K = - C B
def K_ij(height, BC, i, j): 
    hs = height.h_steps
    ws = height.widths
    
    if i == j: #bi : center diag 
        if i == 0 and isinstance(BC, bc.Mixed):

            return -(hs[i+1]**3)/ws[i+1]

        else:
            return -(hs[i]**3)/ws[i] -(hs[i+1]**3)/ws[i+1]
    
    elif i+1 == j: #ci : upper diag
        return (hs[i+1]**3)/ws[i+1]
    
    elif i-1 == j: # ai : lower diag
        return (hs[i]**3)/ws[i]
    
    else:
        return 0


# recursive sequence {Si} 
def get_S(height, BC):
    N = height.N_regions
    
    S = np.zeros(N-2)


    
    off_diag = K_ij(height, BC, N-3, N-2)
    center_diag = K_ij(height, BC, N-2, N-2)
    S[N-3] = -off_diag / center_diag
    
    for k in range(N-4, -1, -1):
        off_diag_succ = off_diag #=K_ij(height, BC, k+1, k+2)
        
        off_diag = K_ij(height, BC, k, k+1)
        center_diag = K_ij(height, BC, k+1, k+1)
        
        S[k] = -off_diag / (center_diag + S[k+1] * off_diag_succ)
        
        
    S_prod = np.zeros(N-1)
    S_prod[0] = S[0]
    
    for k in range(1, N-2):
        S_prod[k] = S_prod[k-1]*S[k]
        
    S_prod[N-2]=1

    return S, S_prod

# recursive sequence {Di} = diags of schur inverse
def get_D(height, BC, S):
    N = height.N_regions
    
    D = np.zeros(N-1)
    

    off_diag = K_ij(height, BC, 0, 1)
    center_diag = K_ij(height, BC, 0, 0)
    D[0] = 1/(center_diag + off_diag * S[0])
    
    
    for k in range(1, N-2):
        off_diag_pred = off_diag
        off_diag = K_ij(height, BC, k, k+1)
        center_diag = K_ij(height, BC, k, k)
        
        D[k] = (1 - off_diag_pred*D[k-1]*S[k-1]) / (center_diag + off_diag * S[k])
        

    off_diag_pred = off_diag
    center_diag = K_ij(height, BC, N-2, N-2)
    D[N-2] = (1 - off_diag_pred*D[N-3]*S[N-3]) / center_diag 
    
    return D

# (P_extrema, P_slopes) -> [p(x)] over domain Nx
def make_ps(height, BC, slopes, extrema):
    ps = np.zeros(height.Nx)
    x0 = height.x0
    k = 0
    x_k = x0
    slope_k = slopes[0]
    
    if isinstance(BC, bc.Fixed):
        p_k = BC.p0
    else:
        p_k = extrema[0]  - slope_k * height.widths[0]
    
    for i in range(height.Nx):
        ps[i] = slope_k*(height.xs[i]-x_k) + p_k
        
        if i >= height.i_peaks[k+1] and k < height.N_regions-1:
            k+= 1
            x_k = height.x_peaks[k]
            p_k = extrema[k-1]
            slope_k = slopes[k]

    return ps

