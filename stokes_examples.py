# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 18:10:34 2025

@author: sarah
"""
import numpy as np
from stokes_heights import PWLinear

class BFS(PWLinear):
    def __init__ (self, args, U, Q, Re, N):
        h_in, h_out, l_in, l_out = args
        x0 = 0
        xf = l_in+l_out
  
        x_peaks = [0, l_in, l_in+l_out]
        y0 = 0
        yf = max(h_in,h_out) 
        y_peaks=[[0,h_in],[h_in,h_out],[h_out,0]]
        namestr= f'BFS_hin{h_in}hout{h_out}lin{l_in}lout{l_out}_U{U}_Q{Q}_Re{Re}'
        super().__init__(x0, xf, y0, yf, N, U, Q, Re,namestr, x_peaks, y_peaks)

        
class BFS_pwl(PWLinear):
    def __init__ (self, args, U, Q, Re, N):
        H, h, L, delta = args
        
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        
        x_peaks = [x0, (L-delta)/2, L/2, (L+delta)/2, xf]
        
        y_peaks=[[0,H],[H,H],[h + (H-h)/2, h + (H-h)/2], [h,h],[h,0]]   
        namestr= f'pwlBFS_H{H}h{h}L{L}d{delta}_U{U}_Q{Q}_Re{Re}'
        super().__init__(x0, xf, y0, yf, N, U, Q, Re,namestr, x_peaks, y_peaks)

class Logistic(PWLinear):
    def __init__(self,args, U, Q, Re, N):
        H, h, L, delta = args
        self.delta = delta
        self.h = h
        self.H = H
        self.center = L//2

        x0 = 0
        xf = L

        x_peaks = [x0 + i/N for i in range (int(1 + L*N))]
        
        y0 = 0
        yf = max(H,h)
        y_peaks_L = [0] + [self.h_fun(x) for x in x_peaks[:-1]]
        y_peaks_R =  [self.h_fun(x) for x in x_peaks[1:]] + [0]

        y_peaks = [[y_L, y_R] for y_L,y_R in np.stack((y_peaks_L,y_peaks_R), axis=1)] 

        
        namestr= f'logistic_H{H}L{L}d{delta}_U{U}_Q{Q}_Re{Re}'
        super().__init__(x0, xf, y0, yf, N, U, Q, Re,namestr, x_peaks, y_peaks)   

    def h_fun(self, x):
        return self.h + (self.H-self.h) / ( 1 + np.exp(self.delta*(x-self.center)))
     
class Sinusoid(PWLinear): #aka sinusoid 2 in reynolds
    def __init__(self, args, U, Q, Re, N):
        H, delta, k,  l, L = args
        
        self.H=H
        self.k = k
        self.delta = delta
        self.l = l

        self.period = np.pi * k / l
        
        x0 = -L
        xf = L
        
        y0 = 0
        yf = H*(1+delta)
        
        x_peaks = [x0 + i/N for i in range (int(1 + 2*L*N))]
        h_peaks = [self.h_fun(x) for x in x_peaks]
        
        y_peaks_L = [0] + h_peaks[:-1]
        y_peaks_R = h_peaks[1:] + [0]

        y_peaks = [[y_L, y_R] for y_L,y_R in np.stack((y_peaks_L,y_peaks_R), axis=1)] 
       
        
        namestr= f'sinusoid_H{H}L{L}d{delta}k{k}_U{U}_Q{Q}_Re{Re}'
        super().__init__(x0, xf, y0, yf, N, U, Q, Re,namestr, x_peaks, y_peaks)


    def h_fun(self,x):
        if abs(x) > self.l:
            if self.k%2 == 0:
                return self.H *(1+self.delta)
            else:
                return self.H *(1-self.delta)
        else:
            return self.H *(1+self.delta*np.cos(self.period * x))
        
          