# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:17:44 2024

@author: sarah
"""

import numpy as np

from scipy.sparse.linalg import LinearOperator 
from scipy.sparse.linalg import gmres
import reyn_boundary as bc




def make_ps(height, BC, cs):

    slopes = height.slopes
    hs = height.hs
    x_peaks = height.x_peaks
    xs = height.xs
    
    ps = np.zeros(height.Nx)
    k = 0
    cq = cs[0] #flux = -cq/12 #/visc
    cu = 6 * BC.U #*visc


    for i in range(height.Nx):
        
        if xs[i] > x_peaks[k+1]:

            k += 1
        
        if slopes[k] != 0:
            dhdx = slopes[k]
            h = hs[i]
            ps[i] = -(cq/2 *h**-2 + cu/h)/dhdx + cs[k+1]
        else: 
            dx = xs[i] - x_peaks[k]
            h = hs[i]
            ps[i] = dx* (cq * h**-3 + cu * h**-2) + cs[k+1]

    return ps


def make_rhs(height, BC):
    N = height.N_regions
    hs = height.h_peaks
    slopes = height.slopes
    widths = height.widths
    
    rhs = np.zeros(N+1)
    c = 6*BC.U #*visc
    
    
    if isinstance(BC, bc.Fixed):
        if slopes[0] != 0:
            rhs[0] = c / (hs[0,1] * slopes[0]) + BC.p0
        else: 
            rhs[0] = BC.p0 
    
    elif isinstance(BC, bc.Mixed):
        rhs[0] = -12 * BC.Q


    for i in range(1, N): # was only until N-1
        
        if slopes[i] != 0 and slopes[i-1] != 0:
            rhs[i] = c *(1/(hs[i,1] * slopes[i]) - 1/(hs[i,0] * slopes[i-1]))
            
        elif slopes[i] != 0 and slopes[i-1] == 0:
            rhs[i] = c * (1/(hs[i,1] * slopes[i]) + hs[i,0]**-2 * widths[i-1])

        elif slopes[i] == 0 and slopes[i-1] != 0:
            rhs[i] = -c /(hs[i,0] * slopes[i-1])

        else: 
            rhs[i] = c * hs[i,0]**-2 * widths[i-1]
            
            
    if slopes[N-1] != 0:
        rhs[N] = c / (hs[N,0]*slopes[N-1]) + BC.pN
    else:
        rhs[N] = -c * widths[N-1]/(hs[N,0]**2) + BC.pN
    
  
    return rhs

#---------------------------------------------------------------------------------
def gmres_solve(height, BC):

    rhs = make_rhs(height, BC)
    linOp = pwlLinOp(height,BC)
    sol_coefs, exit_code = gmres(linOp, rhs, rtol=1e-8)
        
    if exit_code != 0:
        raise Exception('gmres did not converge')

    ps_1D = make_ps(height, BC, sol_coefs)

    return ps_1D

class pwlLinOp(LinearOperator):
    def __init__(self, height,BC):

        self.N = height.N_regions
        self.h_peaks = height.h_peaks
        self.widths = height.widths
        self.slopes = height.slopes
        self.BC = BC
        self.shape = (self.N+1,self.N+1)
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(self.N+1)
    
    def _matvec(self, v):
        hs = self.h_peaks
        slopes = self.slopes
        widths = self.widths
        N = self.N
        mv = np.zeros(N+1)
        cq = v[0]
        
        
        if isinstance(self.BC, bc.Fixed):
            if slopes[0] != 0:
                mv[0] = v[1] - cq/(2 * hs[0,1]**2 * slopes[0])
            else:
                mv[0] = v[1]
                 
        elif isinstance(self.BC, bc.Mixed):
            mv[0] = cq
            
             
        if slopes[N-1] != 0:
            mv[N] = v[N] - cq/(2 * hs[N,0]**2 * slopes[N-1])
        else:
            mv[N] = v[N] + cq  * widths[N-1]/(hs[N,0]**3)
 
 
        for i in range(1, N):
            if slopes[i] != 0 and slopes[i-1] != 0:
                mv[i] = -v[i] + v[i+1] - cq/2 * (1/(hs[i,1]**2 * slopes[i]) - 1/(hs[i,0]**2 * slopes[i-1]))
            
            elif slopes[i] != 0 and slopes[i-1] == 0:
                mv[i] = -v[i] + v[i+1] - cq * (1/(2* hs[i,1]**2 * slopes[i]) + hs[i,0]**-3 * widths[i-1])

            elif slopes[i] == 0 and slopes[i-1] != 0:
                mv[i] = -v[i] + v[i+1] + cq * (1/(2*hs[i,0]**2 * slopes[i-1]))
                
            else:
                mv[i] = -v[i] + v[i+1] - cq * hs[i,0]**-3 * widths[i-1]
           
        return mv

#---------------------------------------------------------------------------------
def schur_solve(height, BC):
    rhs = make_rhs(height, BC)
    dinv_c, dinv_rhs = make_dinv_c_rhs(height,rhs)
    # d_inv_ij = {-1 : i<=j, 0 : i > j}

    N = height.N_regions
    
    cs = np.zeros(N+1) 
    

    cq = rhs[0]  #cq
    
    cs[0] = cq
    
    # cp = (-DinvC + Dinv) * rhs
    for i in range(1,N+1): #cp0 to cpN-1
        cp_i = dinv_c[i-1] * cq + dinv_rhs[i-1]
        cs[i] = cp_i
    
    ps = make_ps(height, BC, cs)
    return ps
    
        

def make_dinv_c_rhs(height, rhs): # make [-D_inv @ C] and [D_inv @ rhs]
    # d_inv_ij = {-1 : i<=j, 0 : i > j}
    slopes = height.slopes  
    widths = height.widths
    hs = height.h_peaks
    N = height.N_regions
                         
    dinv_c = np.zeros(N)  
    dinv_rhs = np.zeros(N) 
    # e: dP bc
    psum_c=0
    psum_rhs=0
    
    if slopes[N-1] != 0:
        psum_c = 1/(2 * hs[N,0]**2 * slopes[N-1])
    else:
        psum_c = -1  * widths[N-1]/(hs[N,0]**3)
    dinv_c[N-1] = psum_c   
    
    psum_rhs = rhs[N] 
    dinv_rhs[N-1] = psum_rhs
    

    for k in range(N-2, -1, -1):
        i=k+1
        if slopes[i] != 0 and slopes[i-1] != 0:
            psum_c +=  - 1/2 * (1/(hs[i,1]**2 * slopes[i]) - 1/(hs[i,0]**2 * slopes[i-1]))
        
        elif slopes[i] != 0 and slopes[i-1] == 0:
            psum_c += - (1/(2* hs[i,1]**2 * slopes[i]) + hs[i,0]**-3 * widths[i-1])

        elif slopes[i] == 0 and slopes[i-1] != 0:
            psum_c += (1/(2*hs[i,0]**2 * slopes[i-1]))
            
        else:
            psum_c +=  - hs[i,0]**-3 * widths[i-1]
        
        dinv_c[k] = psum_c
        
        psum_rhs -= rhs[i]
        dinv_rhs[k] = psum_rhs
    
    return dinv_c, dinv_rhs

 
    
    
