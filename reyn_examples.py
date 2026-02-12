# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

from reyn_heights import PWC_Height,  PWL_Height, LogisticHeight,SinusoidalHeight,  SinusoidalHeight_2
# -----------------------------------------------------------------------------------------------------------------------------------
#   ____________
#  |_____       |
#        |______|
#


class BFS(PWC_Height):
    def __init__(self, args, N):
        h_in, h_out, l_in, l_out = args
        x0 = 0
        xf = l_in+l_out
        N_regions = 2
        x_peaks = np.asarray([0, l_in, xf], float)
        h_peaks = np.asarray([[h_in, h_in], [h_in, h_out], [h_out, h_out]], float)
        namestr = ''#f'BFS_H{int(H)}L{int(xf)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)


# -----------------------------------------------------------------------------------------------------------------------------------
#   ____________
#  |_____       |
#        \______|
#

        
class BFS_pwl(PWL_Height):
    def __init__(self, args, N):
        H, h, L, delta = args

        x0 = 0
        xf = L
        N_regions = 3
        x_peaks = np.asarray([x0, (L-delta)/2, (L+delta)/2, xf], float)
        h_peaks = np.asarray(
            [[0, H], [H, H], [h,h], [h, 0]], float)
        namestr = ''#f'dBFS_H{int(H)}L{int(xf)}_d{int(delta)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)

#----------------------------------------------------------------------------------------
# NOT PIECEWISE LINEAR EXAMPLES...
#-----------------------------------------------------------------------------------------------------------------------------------

class Logistic(LogisticHeight):
    def __init__(self, args, N):
        H_in, H_out, L, delta = args

        # x0 = -L//2
        # xf = L//2
        # center = 0
        x0 = 0
        xf = L
        center = L//2
        namestr = '' #f'Logistic{delta:.2f}H{H:.2f}'
        super().__init__(x0, xf, N, H_in, H_out, center, delta, namestr)

class Sinusoid(SinusoidalHeight): # Reynolds exact solution from Takeuchi & Gu 2019
    def __init__(self, args, N):
        x0 = 0
        
        
        # xf = args[2]
        
        # k = 2*np.pi/xf
        k = args[2] * 2*np.pi
        xf = args[3]

        H = args[0]
        delta = args[1]

        
        self.H = H
        self.delta = delta
        self.k= k
        self.Nx = N * xf
        namestr = ''
        # namestr = f'Sinusoid_H{H}h{h}'
        super().__init__(x0, xf, N, H, delta, k, namestr)

class Sinusoid_2(SinusoidalHeight_2): #dh/dx to 0 at inlet and outlet for compare to Stokes
    def __init__(self, args, N):
        # k = 2*np.pi/xf
        H = args[0]
        delta = args[1]
        k = args[2]
        l = args[3]
        L = args[4]
        
        x0 = -L
        xf = L
    
        
        self.H = H
        self.delta = delta
        self.k= k
        self.Nx = N * xf
        namestr = ''
        # namestr = f'Sinusoid_H{H}h{h}'
        super().__init__(x0, xf, N, H, delta, k, l, L, namestr)        

