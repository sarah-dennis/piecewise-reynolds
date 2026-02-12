#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""
import numpy as np
import reyn_boundary as bc


def fd_solve(height, BC):
    rhs = make_rhs(height, BC)
    mat = make_mat(height, BC)
    ps_1D = np.linalg.solve(mat, rhs)
    return ps_1D

# Reynolds rhs
def make_rhs(height, BC):
    rhs = np.zeros(height.Nx) 
    rhs = 6 * BC.U * height.hxs #*visc
    
    if isinstance(BC, bc.Fixed):
        rhs[0] = BC.p0
        
    elif isinstance(BC, bc.Mixed): 
        h0 = height.hs[0]
        rhs[0] =  -12*BC.Q/h0**3 + 6*BC.U/h0**2  #*visc
        
    rhs[-1] = BC.pN
    return rhs

def make_mat(height, BC):
    N = height.Nx    
    D_lower = np.zeros(N)
    D_center = np.zeros(N)
    D_upper = np.zeros(N)
    
    for i in range(N): #space: xs[i]
        # Find h(t,x) at x = [xs[i-1], xs[i], xs[i+1]] = [hl, hc, hr]
        hl = height.hs[(i-1) % N]
        hc = height.hs[i % N]   
        hr = height.hs[(i+1) % N]

        D_center[i] = -(hr**3 + 2*hc**3 + hl**3) 
        D_upper[i] = (hr**3 + hc**3)
        D_lower[i] = (hl**3 + hc**3)

        
    # combine as upper, middle, lower diagonals
    D = np.diagflat(D_center) + np.diagflat(D_lower[1:N], -1) + np.diagflat(D_upper[0:N-1], 1)
    D /= (2*height.dx**2)
    
    if isinstance(BC, bc.Fixed):
        D[0,0] = 1
        D[0,1] = 0
        
    elif isinstance(BC, bc.Mixed): 
        D[0,0] = -3/(2*height.dx)
        D[0,1] = 4/(2*height.dx)
        D[0,2] = -1/(2*height.dx)
        
    D[N-1,N-1] = 1
    D[N-1,N-2] = 0

    return D
