 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:28:42 2024

@author: sarahdennis
"""
import reyn_pressure_finDiff as fd
import reyn_pressure_pwl as pwl
import reyn_pressure_pwc as pwc

import numpy as np
class Pressure:
    def __init__(self, ps_1D, ps_2D):
        # initializing for adjusted solutions generates both ps_1D and ps_2D
    
        self.ps_1D = ps_1D
        
        self.ps_2D = ps_2D 
    
    def get_dP(self, height):
        if self.ps_2D is None:
            dP =  self.ps_1D[0] - self.ps_1D[-1]

        else:
            ps_2D = np.nan_to_num(self.ps_2D)
            dP = (sum(ps_2D[:,0])/height.hs[0] - sum(ps_2D[:,-1])/height.hs[-1])*height.dy
        return dP

class Reyn_Pressure(Pressure):
    def __init__(self, height, ps_1D):
        ps_2D = self.make_2D_ps(height, ps_1D)
        
        # TODO: if running timing, set ps_2D = None
        # ps_2D = None 
        
        super().__init__(ps_1D, ps_2D)
        
    def make_2D_ps(self, height, ps_1D): # p(x,y) = p(x) 
         ps_2D = np.zeros((height.Ny, height.Nx))
         
         for i in range(height.Nx):
             for j in range(height.Ny):
                 
                 y = height.ys[j]
                 if y <= height.hs[i]:
                     ps_2D[j,i] = ps_1D[i]
                 else:
                     ps_2D[j,i] = None
                
         return ps_2D     
                    
class FinDiff_ReynPressure(Reyn_Pressure):
    def __init__(self, height, BC):
        ps_1D = fd.fd_solve(height, BC)
        
        super().__init__(height, ps_1D)


class PwlGMRes_ReynPressure(Reyn_Pressure):
    def __init__(self, height, BC):
            

            
        ps_1D = pwl.gmres_solve(height, BC)
        super().__init__(height,ps_1D)
        
class PwlSchur_ReynPressure(Reyn_Pressure):
    def __init__(self, height, BC):
            
        ps_1D = pwl.schur_solve(height, BC)
        super().__init__(height,  ps_1D)
    
class PwcSchur_ReynPressure(Reyn_Pressure):
    def __init__(self, height, BC): 
        
        ps_1D = pwc.schur_solve(height, BC)

        super().__init__(height,ps_1D)
        
class PwcSchur_parallel_ReynPressure(Reyn_Pressure):
    def __init__(self, height, BC): 
        
        ps_1D = pwc.schur_solve_parallel(height, BC)


        super().__init__(height, ps_1D)
