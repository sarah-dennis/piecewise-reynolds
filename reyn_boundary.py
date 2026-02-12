# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 09:52:29 2025

@author: sarah
"""
#reynolds

class BoundaryCondition: 
    def __init__(self, U):
        self.pN = 0
        self.U = U
        
         
class Mixed(BoundaryCondition):
    def __init__(self, U, Q):
        self.Q = Q
        super().__init__(U)
        
class Fixed(BoundaryCondition):
    def __init__(self, U, dP):
        self.p0 = dP
        self.dP = dP
        super().__init__(U)

