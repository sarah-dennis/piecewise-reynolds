#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator 
from scipy.sparse.linalg import bicg
import reyn_boundary as bc


def lu_solve(height, BC):
    rhs = make_rhs(height, BC)

    mat = make_mat(height, BC)

    ps_1D = np.linalg.solve(mat, rhs)

    return ps_1D

def bicg_solve(height, BC):
    rhs = make_rhs(height, BC)
    linOp = fdLinOp(height, BC)
    # mat = make_mat(height, BC)
    ps_1D, flag = bicg(linOp, rhs)
    # ps_1D, flag = bicg(mat, rhs)
    if flag:
        raise Exception('biCG did not converge')
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
        
    else: # isinstance(BC, bc.Mixed): 
        D[0,0] = -3/(2*height.dx)
        D[0,1] = 4/(2*height.dx)
        D[0,2] = -1/(2*height.dx)
        
    D[N-1,N-1] = 1
    D[N-1,N-2] = 0

    return D

class fdLinOp(LinearOperator):
    def __init__(self, height, BC):

        self.N = height.Nx
        self.BC = BC
        self.height = height
        self.shape = (self.N,self.N)
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(self.N)
        
    def _matvec(self, v):
       
        mv = np.zeros(self.N)
        
            
        if isinstance(self.BC, bc.Fixed):
            #D[0,0] = 1
            #D[0,1] = 0
            mv[0] = v[0]
            
        else: # isinstance(self.BC, bc.Mixed): 
            # D[0,0] = -3/(2*height.dx)
            # D[0,1] = 4/(2*height.dx)
            # D[0,2] = -1/(2*height.dx)
            mv[0] = (-3*v[0] + 4*v[1] -1*v[2])/(2*self.height.dx)
         
        mv[-1] = v[-1]
        
        
        for i in range(1,self.N-1): #space: xs[i]
            # Find h(t,x) at x = [xs[i-1], xs[i], xs[i+1]] = [hl, hc, hr]
            hl = self.height.hs[(i-1)]
            hc = self.height.hs[i]   
            hr = self.height.hs[(i+1)]
    
            D_center = -(hr**3 + 2*hc**3 + hl**3)   # row i, col j   
            D_upper = (hr**3 + hc**3)               #        col j+1 
            D_lower = (hl**3 + hc**3)               #        col j-1
            
            mv[i] = (D_center*v[i] + D_upper*v[i+1] + D_lower*v[i-1])/(2*self.height.dx**2)

           
        return mv

    def _rmatvec(self, v):
       
        rmv = np.zeros(self.N)
        
            
        if isinstance(self.BC, bc.Fixed):
            bc_0 = 1 
            bc_1 = 0
            bc_2 = 0
            
        else: # isinstance(self.BC, bc.Mixed): 

            bc_0 = -3/(2*self.height.dx)
            bc_1 = 4/(2*self.height.dx)
            bc_2 = -1/(2*self.height.dx)
        
        #i=0
        hl_1 = self.height.hs[0]
        hc_1 = self.height.hs[1]   
        D_lower_1 = (hl_1**3 + hc_1**3)/(2*self.height.dx**2)
        rmv[0] = bc_0*v[0] + D_lower_1*v[1]
        
        #i=1
        hl_1  = self.height.hs[0]
        hc_1 = self.height.hs[1]
        hr_1 = self.height.hs[2]

        hl_2 = hc_1
        hc_2 = hr_1
        
        D_center_1 = -(hr_1**3 + 2*hc_1**3 + hl_1**3)/(2*self.height.dx**2)
        D_lower_2 = (hl_2**3 + hc_2**3)/(2*self.height.dx**2)
        
        rmv[1] = bc_1*v[0] + D_center_1*v[1] + D_lower_2*v[2]         
            
        #i=2
        hr_2 = self.height.hs[3]
    
        hl_3 = hc_2
        hc_3 = hr_2

        D_center_2 = -(hr_2**3 + 2*hc_2**3 + hl_2**3)/(2*self.height.dx**2)   # row i, col j   
        D_upper_1 = (hr_1**3 + hc_1**3) /(2*self.height.dx**2)              #        col j+1 
        D_lower_3 = (hl_3**3 + hc_3**3)/(2*self.height.dx**2)               #        col j-1
        
        rmv[2] = bc_2*v[0] + D_upper_1*v[1] + D_center_2*v[2] + D_lower_3*v[3]
        
        
        #i=N-1
        hr_ = self.height.hs[-1]
        hc_  = self.height.hs[-2]
        
        D_upper_ = (hr_**3 + hc_**3)/(2*self.height.dx**2)
        
        rmv[-1] = v[-1] + D_upper_*v[-2]
        
        
        for i in range(3,self.N-1): #space: xs[i]
            # Find h(t,x) at x = [xs[i-1], xs[i], xs[i+1]] = [hl, hc, hr]
            hl = self.height.hs[(i-1)]
            hc = self.height.hs[i]   
            hr = self.height.hs[(i+1)]
    
            hr_minus = hc
            hc_minus = hl
            hl_plus = hc
            hc_plus = hr
            
            D_center = -(hr**3 + 2*hc**3 + hl**3)   # row i, col j   
            D_upper_minus = (hr_minus**3 + hc_minus**3)               #        col j+1 
            D_lower_plus = (hl_plus**3 + hc_plus**3)               #        col j-1
            
            rmv[i] = (D_upper_minus*v[i-1] + D_center*v[i] + D_lower_plus*v[i+1])/(2*self.height.dx**2)

           
        return rmv
