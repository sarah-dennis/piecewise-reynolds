 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 08:18:07 2023

@author: sarahdennis
"""
import numpy as np
from domain import Height
# import graphics
#------------------------------------------------------------------------------
# PWL Height
#------------------------------------------------------------------------------
class PWL_Height(Height):
    def __init__(self, x0, xf, N, N_regions, x_peaks, h_peaks, filestr):
        self.h_peaks = h_peaks
        self.x_peaks = x_peaks
        self.N_regions = N_regions # = len(h_peaks)-1
        
        hs, self.slopes, self.widths, i_peaks = self.make_hs(x0, xf, N, x_peaks, h_peaks)
        
        y0 = 0
        yf = max(hs)  
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)
        
    def make_hs(self, x0, xf, N, x_peaks, h_peaks):
        slopes = np.zeros(self.N_regions)
        widths = np.zeros(self.N_regions)
        i_peaks = np.zeros(self.N_regions+1, dtype=int)
        for r in range(self.N_regions):            
            slopes[r] = (h_peaks[r+1,0] - h_peaks[r,1])/(x_peaks[r+1] - x_peaks[r])
            
        Nx = int((xf-x0)*N + 1)
        hs = np.zeros(Nx)
        dx = 1/N
        
        r = 0
        for i in range(Nx):
            xi = x0 + i*dx
            
            if xi > x_peaks[r+1] and r+1 < self.N_regions:
                r +=1
            
            i_peaks[r+1] = i    
            widths[r] = xi - x_peaks[r]
            hs[i] = h_peaks[r,1] + slopes[r] * (xi - x_peaks[r])
        # print(hs)
        
        return  hs, slopes, widths, i_peaks
    
def make_PWL(height):  
    
    N = height.N
    N_regions = height.Nx-1
    x_peaks = height.xs
    h_peaks = np.zeros((height.Nx, 2))

    h_peaks[0,0] = 0
    
    
    for i in range(height.Nx-1):
        
        h_peaks[i,1] =  height.hs[i]
        h_peaks[i+1,0] =  height.hs[i+1]
        
    
    h_peaks[height.Nx-1,1] = 0

    new_height = PWL_Height(height.x0, height.xf, N, N_regions, x_peaks, h_peaks, '')
    
    return new_height
    
    
#------------------------------------------------------------------------------
# PWC Height
#------------------------------------------------------------------------------
class PWC_Height(Height):#(PWL_Height):
    def __init__(self, x0, xf, N, N_regions, x_peaks, h_peaks, filestr):

        #solver requires minimum 3 regions
        while N_regions <3 :
            x_peak_mid = x_peaks[-1] - x_peaks[-2]/2
            x_peak_end = x_peaks[-1]
            x_peaks[-1] = x_peak_mid
            x_peaks = np.append(x_peaks, x_peak_end)
            h_peaks = np.append(h_peaks, h_peaks[-1])
            h_peaks = np.reshape(h_peaks, (N_regions+2,2))
            N_regions +=1 
        self.N_regions = N_regions
        self.x_peaks=x_peaks
        self.h_peaks=h_peaks
        self.h_steps = h_peaks[:,1][:-1] 
        self.slopes =np.zeros(self.N_regions)
        
        hs, self.widths, i_peaks = self.make_hs(x0, xf, N, x_peaks, h_peaks)
        y0 = 0
        yf = max(hs)  
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)
        
    def make_hs(self, x0, xf, N, x_peaks, h_peaks):
        widths = np.zeros(self.N_regions)
        i_peaks = np.zeros(self.N_regions+1, dtype=int)
       
        Nx = int((xf-x0)*N + 1)
        hs = np.zeros(Nx)
        dx = 1/N
        
        r = 0
        for i in range(Nx):
            xi = x0 + i*dx
            if xi > x_peaks[r+1] and r+1 < self.N_regions:
                r +=1
            i_peaks[r+1] = i    
            widths[r] = xi - x_peaks[r]
            hs[i] = h_peaks[r,1]
        return  hs, widths, i_peaks

def make_PWC(height):

    N = height.N
    N_regions = height.Nx-1
    x_peaks = height.xs
    h_peaks = np.zeros((height.Nx, 2))

    h_peaks[0,0] = 0
    
    
    for i in range(height.Nx-1):
        h_step = (height.hs[i+1] + height.hs[i])/2
        
        h_peaks[i,1] = h_step
        h_peaks[i+1,0] = h_step
        
    
    h_peaks[height.Nx-1,1] = 0
    
    new_height = PWC_Height(height.x0, height.xf, N, N_regions, x_peaks, h_peaks, '')
    # print(new_height.h_steps, height.hs)
    return new_height
 
#------------------------------------------------------------------------------
# Other Height Functions
#------------------------------------------------------------------------------
class SinusoidalHeight(Height):  #exact solution
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, x0, xf, N, H, delta, k, filestr):
        
        self.H = H
        self.k = k
        self.delta=delta
        
       
        y0 = 0
        yf = H *(1+delta)
        
        Nx = (xf-x0)*N + 1
        dx = 1/N
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        
        i_peaks = [0, Nx-1]

        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)

    
    def h_fun(self,x):
        return self.H *(1+self.delta*np.cos(self.k*x))
    
    def h_recip_deriv_fun(self,x):
        num = self.delta * self.k  * np.sin(self.k * x)
        dnm = -self.H * (1 + self.delta * np.cos(self.k*x))**2
        return  num/dnm
    
class SinusoidalHeight_2(Height): # inlet and outlet dh/dx=0
    def __init__(self, x0, xf, N, H, delta, k, l, L, filestr):
        
        self.H = H
        self.k = k
        self.delta = delta
        self.l = l
        self.L = L
        self.period = k*np.pi/l
        
        y0 = 0
        yf = H *(1+delta)
        
        
        Nx = 2*L*N + 1
        dx = 1/N
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        
        i_peaks = [0, Nx-1]
  
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)
  
    
    def h_fun(self,x):
        if abs(x) > self.l:
            if self.k%2 == 0:
                return self.H *(1+self.delta)
            else:
                return self.H *(1-self.delta)
        else:
            return self.H *(1+self.delta*np.cos(self.period * x))
        

class BumpHeight(Height): 
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, x0, xf, N, lam, H, h0, filestr):
        Nx = (xf-x0)*N + 1

        y0 = 0
        
        self.H=H         
        self.x_scale = (xf-x0)/2  
        self.h0 = h0
        self.lam=lam
        
        dx = 1/N
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        yf = np.max(hs)
        i_peaks = [0, Nx-1]
       
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)
  
    def h_fun(self, x):
        return self.H*(1-(self.lam/2)*(1+np.cos(np.pi*x/self.x_scale)))-(self.H-self.h0)
   
class LogisticHeight(Height):
    def __init__(self, x0, xf, N, H, h, center, delta, filestr):
       
        Nx = (xf-x0)*N + 1
        dx = 1/N
        self.h = h
        self.H = H
        self.delta = delta #slope = delta*(H-h)/4
        self.center = center
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])  
        y0 = 0
        yf = np.max(hs)
        i_peaks = [0, Nx-1]
       
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)

    def h_fun(self, x):
        
        return (self.h + (self.H-self.h) / ( 1 + np.exp(self.delta*(x-self.center))))  

#------------------------------------------------------------------------------    
class CircleHeight(Height):
    
    def __init__(self, x0, xf, N, r, dxdr, h0, filestr):
        Nx = int((xf-x0)*N + 1)
        dx = 1/N
        
        self.h0 = h0 #clearance
        self.r = r 
        self.l = ((xf-x0)-2*r)/2
            
        self.dxdr =  dxdr
        
            
        dx = (xf - x0)/(Nx-1)

        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        y0 = 0
        yf = max(hs)
        # print(self.l, self.l*N)
        if self.l != 0:
            i_peaks = np.asarray([0, (self.l+self.dxdr)*N, Nx-1-(self.l+self.dxdr)*N, Nx-1], int)
        else:
            i_peaks = [0, Nx-1]
       
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)

    def h_fun(self, x):
        # return self.h0 + np.sqrt(self.r**2 - (x-self.r)**2)
        if x <-self.r+self.dxdr or x >self.r-self.dxdr:
            return self.h0+self.r - np.sqrt(self.r**2 - (self.r-self.dxdr)**2)
        else:
            return self.h0 + self.r - np.sqrt(self.r**2 - x**2)
        
    
    