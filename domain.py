#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 17:33:35 2023

@author: sarahdennis
"""
import numpy as np
import graphics
#------------------------------------------------------------------------------
# Domain: 2D grid on [x0, xf] and [y0, yf] on Nx and Ny grid points

class Domain:

    def __init__(self, x0, xf, y0, yf, N, dirstr=None, namestr=None):
        self.x0 = x0            # left boundary x = x0
        self.xf = xf            # right boundary x = xf 
        self.y0 = y0            # lower surface y = y0
        self.yf = yf            # upper boundary y = max h = yf-y0
        
        self.N = N  # scale: number of grid points in [0,1]
        self.Nx = int((xf-x0)*N + 1) #total number of x grid points
        self.Ny = int((yf-y0)*N + 1) #total number of y grid points
        
        self.dx = 1/N
        self.dy = 1/N
        
        self.xs = np.linspace(x0, xf, self.Nx)
        self.ys = np.linspace(y0, yf, self.Ny)
        self.dirstr = dirstr
        self.filestr= f"{dirstr}/{namestr}_N{N}"
        

#------------------------------------------------------------------------------
# Domain for Reynolds solver
class Height(Domain):

    def __init__(self, x0, xf, y0, yf, N, hs, i_peaks, namestr):
        dirstr = f"./reyn_examples/{namestr}"
        super().__init__(x0, xf, y0, yf, N, dirstr, namestr) # -> {dx, dy, xs, ys}
        
        self.hs = hs
        self.H_in=hs[0]
        self.H_out=hs[-1]
        self.i_peaks = i_peaks
    
                
        # graphics.plot_2D_multi([self.hxs,self.h2xs,self.h3xs], self.xs, 'Height gradients', ['$h_x$','$h_{xx}$','$h_{xxx}$'], ['$x$','$h_{*}$'])

        self.h_min = min(self.hs)
      
    
# Domain for Stokes solver
class Space(Domain):
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, namestr):
        dirstr = f"./stokes_examples/{namestr}"
        super().__init__(x0, xf, y0, yf, N, dirstr, namestr)
        self.U = U    # velocity at flat boundary 
        # self.visc = 1  # dynamic viscosity (Ns/m^2)
        # self.dens = 1  # density (kg/m^3)

        self.p_ambient = 0 #pressure (outlet) Pa   
        self.flux=Q
        self.Re = Re #

#------------------------------------------------------------------------------
def center_diff(fs, Nx, dx):
    D_lower = -1*np.ones(Nx)
    D_upper = np.ones(Nx)
    D = np.diagflat(D_lower[1:Nx], -1) + np.diagflat(D_upper[0:Nx-1], 1)
    D = D/(2*dx)
    
    fs_dx = D@fs 
    
    # fs_dx[0] = fs_dx[1]
    fs_dx[0]= right_first(dx, fs[0 : 3]) 
    
    # fs_dx[-1] = fs_dx[-2]
    fs_dx[Nx-1]= left_first(dx, fs[Nx-3 : Nx])
    return np.asarray(fs_dx)
      
def center_second_diff(fs, Nx, dx):
    D_lower = np.ones(Nx-1)
    D_upper = np.ones(Nx-1)
    D_center = -2*np.ones(Nx)
    D = np.diagflat(D_lower, -1) +  np.diagflat(D_center, 0)+ np.diagflat(D_upper, 1)

    D = D/(dx**2)
    fs_dxx = D@fs 
    
    fs_dxx[0]= right_second(dx, fs[0 : 4]) #(2*fs[0]-5*fs[1]+4*fs[2]-fs[3])/(dx**2)
    fs_dxx[Nx-1]= left_second(dx, fs[Nx-4 : Nx]) #(2*fs[Nx-1]-5*fs[Nx-2]+4*fs[Nx-3]-fs[Nx-4])/(dx**2)
    
    return np.asarray(fs_dxx)

 # return np.dot([-1, 2, 0, -2, 1], u2W2E)/(2*dx**3)
def center_third_diff(fs, Nx, dx):
    # D_lllower = 1*np.ones(Nx-3)
    # D_llower = -8*np.ones(Nx-2)
    # D_lower = 13*np.ones(Nx-1)
    # D_upper = -13*np.ones(Nx-1)
    # D_uupper = 8*np.ones(Nx-2)
    # D_uuupper = -1*np.ones(Nx-3)

    # D = np.diagflat(D_lllower, -3) +np.diagflat(D_llower, -2) +np.diagflat(D_lower, -1) + np.diagflat(D_upper, 1) +np.diagflat(D_uupper, 2)+np.diagflat(D_uuupper, 3)

    # D = D/(8*dx**3)
    # fs_dxxx = D@fs 
    
    D_llower = -1*np.ones(Nx-2)
    D_lower = 2*np.ones(Nx-1)
    D_upper = -2*np.ones(Nx-1)
    D_uupper = 1*np.ones(Nx-2)
    # D_center = -2*np.ones(Nx)
    D = np.diagflat(D_llower, -2) +np.diagflat(D_lower, -1) + np.diagflat(D_upper, 1) +np.diagflat(D_uupper, 2)

    D = D/(2*dx**3)
    fs_dxxx = D@fs 
    
    fs_dxxx[0]= right_third(dx, fs[0 : 5]) #(2*fs[0]-5*fs[1]+4*fs[2]-fs[3])/(dx**2)
    fs_dxxx[1]= right_third(dx, fs[1 : 6]) #(2*fs[0]-5*fs[1]+4*fs[2]-fs[3])/(dx**2)
    # fs_dxxx[2]= right_third(dx, fs[2 : 7]) #(2*fs[0]-5*fs[1]+4*fs[2]-fs[3])/(dx**2)
    
    fs_dxxx[Nx-1]= left_third(dx, fs[Nx-5 : Nx]) #(2*fs[Nx-1]-5*fs[Nx-2]+4*fs[Nx-3]-fs[Nx-4])/(dx**2)
    fs_dxxx[Nx-2]= left_third(dx, fs[Nx-6 : Nx-1]) #(2*fs[Nx-1]-5*fs[Nx-2]+4*fs[Nx-3]-fs[Nx-4])/(dx**2)
    # fs_dxxx[Nx-3]= left_third(dx, fs[Nx-7 : Nx-2]) #(2*fs[Nx-1]-5*fs[Nx-2]+4*fs[Nx-3]-fs[Nx-4])/(dx**2) 
    
    return np.asarray(fs_dxxx)
# [1, -4, 6, -4, 1], u2W2E)/(dx**4)
def center_fourth_diff(fs, Nx, dx):
    D_llower = 1*np.ones(Nx-2)
    D_lower = -4*np.ones(Nx-1)
    D_center = 6*np.ones(Nx)
    D_upper = -4*np.ones(Nx-1)
    D_uupper = 1*np.ones(Nx-2)
    
    D = np.diagflat(D_llower, -2) +np.diagflat(D_lower, -1) + np.diagflat(D_center, 0)+np.diagflat(D_upper, 1) +np.diagflat(D_uupper, 2)

    D = D/(dx**4)
    fs_dxxxx = D@fs 
    
    fs_dxxxx[0]= right_fourth(dx, fs[0 : 6]) #(2*fs[0]-5*fs[1]+4*fs[2]-fs[3])/(dx**2)
    fs_dxxxx[1]= right_fourth(dx, fs[1 : 7]) #(2*fs[0]-5*fs[1]+4*fs[2]-fs[3])/(dx**2)
    fs_dxxxx[Nx-1]= left_fourth(dx, fs[Nx-6 : Nx]) #(2*fs[Nx-1]-5*fs[Nx-2]+4*fs[Nx-3]-fs[Nx-4])/(dx**2)
    fs_dxxxx[Nx-2]= left_fourth(dx, fs[Nx-7 : Nx-1]) #(2*fs[Nx-1]-5*fs[Nx-2]+4*fs[Nx-3]-fs[Nx-4])/(dx**2)
    
    return np.asarray(fs_dxxxx)

#----------------------------------------------------------------------------

def right_first_O1(dx, uuE):
    return  np.dot([-1, 1], uuE)/dx

def right_first(dx, uu2E):
    # return (-3*u + 4*uE - u2E)/(2*dx)
    return np.dot([-3, 4, -1], uu2E) /(2*dx)

def right_first_O4(dx, uu4E):
    # return (-3*u + 4*uE - u2E)/(2*dx)
    return np.dot([-25,48,-36,16,-3], uu4E)/(12*dx)
    
def right_second_O1(dx, uu2E):
    return np.dot([1, -2, 1], uu2E)/(dx**2)

def right_second(dx, uu3E):
    # return (2*u - 5*uE + 4*u2E - u3E)/(dx**2)
    return np.dot([2, -5, 4, -1], uu3E)/(dx**2)

def right_second_O4(dx, uu5E):
    # return (2*u - 5*uE + 4*u2E - u3E)/(dx**2)
    return np.dot([45, -154, 214, -156, 61,-10], uu5E)/(12*dx**2)

def right_third(dx, uu4E):
    # return (-5*u + 18*uE -24*u2E + 14*u3E -3*u4E)/(2*dx**3)
    return np.dot([-5, 18, -24, 14, -3], uu4E)/(2*dx**3)

def right_third_O4(dx, uu6E):
    # return (-5*u + 18*uE -24*u2E + 14*u3E -3*u4E)/(2*dx**3)
    return np.dot([-49,232,-461,496,-307,104,-15], uu6E)/(8*dx**3)
  
def right_fourth(dx, uu5E):
    return np.dot([3, -14, 26, -24, 11, -2], uu5E)/(dx**4)
  
def right_fourth_O4(dx, uu7W):
    return np.dot([28/3, -111/2, 142, -1219/6, 176, -185/2, 82/3, -7/2], uu7W)/(dx**4)

def right_fifth(dx, uu6E):
    return np.dot([-7, 40, -95, 120, -85, 32, -5], uu6E)/(2*dx**5)    

#----------------------------------------------------------------------------
def left_first_O1(dx, uuW):  
    return np.dot([-1,1], uuW)/dx
    # return (uW-u) /(dx)

def left_first(dx, uu2W):
    # return (3*u - 4*uW + u2W)/(2*dx)    
    return np.dot([1, -4, 3], uu2W) /(2*dx)

def left_first_O4(dx, uu4W):
    # return (-3*u + 4*uE - u2E)/(2*dx)
    return np.dot([3, -16, 36, -48, 25], uu4W)/(12*dx)

def left_second_O1(dx, uu2W):
    # return (2*u - 5*uW + 4*u2W - u3W)/(dx**2)
    return np.dot([1, -2, 1], uu2W)/(dx**2) 

def left_second(dx, uu3W):
    # return (2*u - 5*uW + 4*u2W - u3W)/(dx**2)
    return np.dot([-1, 4, -5, 2], uu3W)/(dx**2)

def left_second_O4(dx, uu5W):
    # return (2*u - 5*uE + 4*u2E - u3E)/(dx**2)
    return np.dot([-10, 61, -156, 214, -154, 45], uu5W)/(12*dx**2)

def left_third(dx, uu4W):
    # return (5*u - 18*uW +24*u2W - 14*u3W +3*u4W)/(2*dx**3)
    return np.dot([3, -14, 24, -18, 5], uu4W)/(2*dx**3)

def left_third_O4(dx, uu6W):
    # return (-5*u + 18*uE -24*u2E + 14*u3E -3*u4E)/(2*dx**3)
    return np.dot([15, -104, 307, -496, 461, -232, 49], uu6W)/(8*dx**3)

def left_fourth(dx, uu5W):
    return np.dot([-2, 11, -24, 26, -14, 3], uu5W)/(dx**4)

def left_fourth_O4(dx, uu7W):
    return np.dot([-7/2, 82/3, -185/2, 176, -1219/6, 142, -111/2, 28/3], uu7W)/(dx**4)

def left_fifth(dx, uu6W):
    return np.dot([5, -32, 85, -120, 95, -40, 7], uu6W)/(2*dx**5)
#------------------------------------------------------------------------------
def center_first(dx, uWE):
    # return (uW - uE)/(2*dx)
    return np.dot([-1, 0, 1], uWE)/(2*dx)

def center_first_O4(dx, u2W2E):
    # return (uW - uE)/(2*dx)
    return np.dot([1, -8, 0, 8, -1], u2W2E)/(12*dx)

def center_second(dx, uWE):
    return np.dot([1,-2, 1],uWE)/(dx**2)
    # return (uE -2*u + uW)/(dx**2)
    
def center_second_O4(dx, u2W2E):
    # return (uW - uE)/(2*dx)
    return np.dot([-1, 16, -30, 16, -1], u2W2E)/(12*dx**2)

def center_third(dx, u2W2E):
    # return (u2E - 2*uE +2*uW -u2W)/(2*dx**3)
    return np.dot([-1, 2, 0, -2, 1], u2W2E)/(2*dx**3)

def center_third_O4(dx, u3W3E):
    # return (uW - uE)/(2*dx)
    return np.dot([1, -8, 13, 0, -13, 8, -1], u3W3E)/(8*dx**3)

def center_fourth(dx, u2W2E):
    return np.dot([1, -4, 6, -4, 1], u2W2E)/(dx**4)

def center_fifth(dx, u3W3E):
    return np.dot([-1, 4, -5, 0, 5,-4, 1], u3W3E)/(2*dx**5)

#------------------------------------------------------------------------------

        
def avg_x(fs):
    # i=1
    # f = (fs[i+1] + fs[i-1])/2 
    # return f   
    return avg_2x(fs)
                   
def avg_2x(fs):
    #-2,-1,0,1,2
    
    i=2
    f = (fs[i+2] + fs[i-2])/2 
    f_E = (f + fs[i+2])/2
    f_W = (f + fs[i-2])/2

    return np.array([f_W, f, f_E], dtype='f8')   

def avg_3x(fs):
    #-3,-2,-1,0,1,2,3
    i=3
    f = (fs[i+3] + fs[i-3])/2 
    f_E = (f + fs[i+3])/2
    f_W = (f + fs[i-3])/2
    f_2E = (f_E + fs[i+3])/2
    f_2W = (f_W + fs[i-3])/2
    return np.array([f_2W, f_W, f, f_E, f_2E]  ) 

def avg_4x(fs):
    #-4,-3,-2,-1,0,1,2,3,4
    i=4
    f = (fs[i+4] + fs[i-4])/2 
    f_E = (f + fs[i+4])/2
    f_W = (f + fs[i-4])/2
    f_2E = (f_E + fs[i+4])/2
    f_2W = (f_W + fs[i-4])/2
    f_3E = (f_2E + fs[i+4])/2
    f_3W = (f_2W + fs[i-4])/2
    return np.array([f_3W, f_2W, f_W, f, f_E, f_2E, f_3E]   )

def avg_5x(fs):
    #-5,-4,-3,-2,-1,0,1,2,3,4,5
    i=5
    f = (fs[i+5] + fs[i-5])/2 
    f_E = (f + fs[i+5])/2
    f_W = (f + fs[i-5])/2
    f_2E = (f_E + fs[i+5])/2
    f_2W = (f_W + fs[i-5])/2
    f_3E = (f_2E + fs[i+5])/2
    f_3W = (f_2W + fs[i-5])/2
    f_4E = (f_3E + fs[i+5])/2
    f_4W = (f_3W + fs[i-5])/2
    return np.array([f_4W,f_3W, f_2W, f_W, f, f_E, f_2E, f_3E,f_4E]  ) 


def avg_6x(fs):
    #-5,-4,-3,-2,-1,0,1,2,3,4,5
    i=6
    f = (fs[i+6] + fs[i-6])/2 
    f_E = (f + fs[i+6])/2
    f_W = (f + fs[i-6])/2
    f_2E = (f_E + fs[i+6])/2
    f_2W = (f_W + fs[i-6])/2
    f_3E = (f_2E + fs[i+6])/2
    f_3W = (f_2W + fs[i-6])/2
    f_4E = (f_3E + fs[i+6])/2
    f_4W = (f_3W + fs[i-6])/2
    f_5E = (f_4E + fs[i+6])/2
    f_5W = (f_4W + fs[i-6])/2
    return np.array([f_5W,f_4W,f_3W, f_2W, f_W, f, f_E, f_2E, f_3E,f_4E,f_5E]  ) 

def avg_7x(fs):
    #-5,-4,-3,-2,-1,0,1,2,3,4,5
    i=7
    f = (fs[i+7] + fs[i-7])/2 
    f_E = (f + fs[i+7])/2
    f_W = (f + fs[i-7])/2
    f_2E = (f_E + fs[i+7])/2
    f_2W = (f_W + fs[i-7])/2
    f_3E = (f_2E + fs[i+7])/2
    f_3W = (f_2W + fs[i-7])/2
    f_4E = (f_3E + fs[i+7])/2
    f_4W = (f_3W + fs[i-7])/2
    f_5E = (f_4E + fs[i+7])/2
    f_5W = (f_4W + fs[i-7])/2
    f_6E = (f_5E + fs[i+7])/2
    f_6W = (f_5W + fs[i-7])/2
    return np.array([f_6W,f_5W,f_4W,f_3W, f_2W, f_W, f, f_E, f_2E, f_3E,f_4E,f_5E,f_6E]) 

def avg_8x(fs):
    #-5,-4,-3,-2,-1,0,1,2,3,4,5
    i=8
    f = (fs[i+8] + fs[i-8])/2 
    f_E = (f + fs[i+8])/2
    f_W = (f + fs[i-8])/2
    f_2E = (f_E + fs[i+8])/2
    f_2W = (f_W + fs[i-8])/2
    f_3E = (f_2E + fs[i+8])/2
    f_3W = (f_2W + fs[i-8])/2
    f_4E = (f_3E + fs[i+8])/2
    f_4W = (f_3W + fs[i-8])/2
    f_5E = (f_4E + fs[i+8])/2
    f_5W = (f_4W + fs[i-8])/2
    f_6E = (f_5E + fs[i+8])/2
    f_6W = (f_5W + fs[i-8])/2
    f_7E = (f_6E + fs[i+8])/2
    f_7W = (f_6W + fs[i-8])/2
    return np.array([f_7W,f_6W,f_5W,f_4W,f_3W, f_2W, f_W, f, f_E, f_2E, f_3E,f_4E,f_5E,f_6E, f_7E]) 