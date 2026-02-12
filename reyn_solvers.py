# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 10:10:36 2026

@author: sarah
"""
import numpy as np
from time import time

import domain as dm
import reyn_boundary as bc
import reyn_velocity as rv
import reyn_pressure as rp
from reyn_control import Reyn_Solution
from reyn_heights import PWC_Height, PWL_Height, make_PWC, make_PWL


class Reynolds_Solver: 
    def __init__(self, Example, BC, args=None):
        self.Example = Example #initialize height = Example(args) in the solver
        self.args = args
        
        self.BC = BC        
        
#----------------------------------------------------------------------------------
    def fd_solve(self, N):
        solver_title = "Reynolds"
        
        height = self.Example(self.args, N)

        t0 = time()
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        pressure = rp.FinDiff_ReynPressure(height, self.BC)
        tf = time()
            
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        t = tf-t0
        
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
#----------------------------------------------------------------------------------
#---------------------------------Piecewise----------------------------------------
#----------------------------------------------------------------------------------
    def pwc_solve(self, N):
        solver_title = "Reynolds"
        
        height = self.Example(self.args, N)

        if not isinstance(height, PWC_Height):
            height = make_PWC(height)
    
        
        t0 = time()
        pressure = rp.PwcSchur_ReynPressure(height, self.BC)
        tf = time()
            
        height.hxs = np.zeros(height.Nx)
        velocity = rv.Reyn_Velocity(height, self.BC,pressure)
        
        t=tf-t0
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
#----------------------------------------------------------------------------------   
   
    def pwc_parallel_solve(self, N):
        solver_title = "Reynolds"
        
        height = self.Example(self.args, N)

        if not isinstance(height, PWC_Height):
            height = make_PWC(height)
        
        t0 = time()
        pressure = rp.PwcSchur_parallel_ReynPressure(height, self.BC)
        tf = time()
        
        velocity = rv.Reyn_Velocity(height, self.BC,pressure)
        

        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        
        t = tf-t0
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution

#----------------------------------------------------------------------------------   
    def pwl_solve(self, N):
        solver_title = "Reynolds" #" Piecewise Linear"
        
        height = self.Example(self.args,N)
        
        if not (isinstance(self.BC, bc.Mixed)): #TODO
            raise TypeError('Only prescribed flux for pwl schur solver')
        
        if not (isinstance(height, PWL_Height) or isinstance(height, PWC_Height)): # PWC is a PWL
            height = make_PWL(height)

        t0 = time()
        pressure = rp.PwlSchur_ReynPressure(height, self.BC)
        tf = time()
        t = tf-t0
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        
        
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
        return pressure, velocity, tf-t0
    
#----------------------------------------------------------------------------------       
    def pwl_gmres_solve(self, N):
        solver_title = "Reynolds" #" Piecewise Linear"
        
        height = self.Example(self.args,N)
        
        t0 = time()
       
        if not (isinstance(self.BC, bc.Mixed)): #TODO
            raise TypeError('Only prescribed flux for pwl gmres solver')
 
       
        if not (isinstance(height, PWL_Height) or isinstance(height, PWC_Height)) :
            height = make_PWL(height)
       
        pressure = rp.PwlGMRes_ReynPressure(height, self.BC)
        tf = time()
        
        # print('pwl gmres time: ', tf-t0)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
            
        t = tf-t0
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
    
    
        
#----------------------------------------------------------------------------------
    def exact_sinusoid_sol(self, N):
        height = self.Example(self.args, N)
       
        ps = np.zeros(height.Nx)
        t0 = time()
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        h_recip_dx = [height.h_recip_deriv_fun(x) for x in height.xs]
        for i in range(height.Nx):
            h = height.hs[i] 
            ps[i] = -6*self.BC.U * (h + height.H)/((height.k * height.H)**2 * (2+height.delta**2)) * h_recip_dx[i]
        
        pressure = rp.Reyn_Pressure(height, ps)
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        tf = time()
        t = tf-t0
        solver_title = "exact solution"
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)
        return solution