# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:15:15 2024

@author: sarah
"""
import numpy as np
import math
from domain import Space

#------------------------------------------------------------------------------
class PWLinear(Space):
    # peak_xs = [x0, ..., xi,..., xf] : x0 < xi < xf
    # peak_ys = [(yf,h_in),...,(hi_left, hi_right),...,(h_out,yf)]
    
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, namestr, x_peaks, y_peaks):
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, namestr)
        # peaks must fall on the grid
        self.x_peaks = x_peaks

        self.y_peaks = y_peaks
        self.N_regions = len(self.x_peaks)-1
        self.make_space()

        self.H_in = y_peaks[0][1] - y0
        self.H_out = y_peaks[-1][0] - y0
        
        self.spacestr = "$Q=%.2f$, $U=%.1f$"%(Q,U)  # for plot title
        
        if self.H_in == 0: # closed cavity --> Q=0, dp=0
            self.dp_in = 0
        else:              # gap entry --> dp ~ Q
            self.dp_in = (self.flux - 0.5*self.U*self.H_in) * (-12 / self.H_in**3)
    
    # 0 : boundary, -1: exterior, 1: interior
    def make_space(self):

        slopes = np.zeros(self.N_regions)
        
        for k in range(self.N_regions):
            dh = self.y_peaks[k+1][0] - self.y_peaks[k][1]
            dx = self.x_peaks[k+1] - self.x_peaks[k]
            slopes[k] = dh/dx
        hs = np.zeros((self.Nx,2))
        grid = np.zeros((self.Ny, self.Nx))
    
        reg = 0        
        i_ref = 0
        for i in range(self.Nx):
            if math.isclose(self.xs[i], self.x_peaks[reg]):
                h_left = self.y_peaks[reg][0]
                h_right = self.y_peaks[reg][1] 
                i_ref = i
                reg +=1
                hs[i] = [h_left,h_right]
            else:
                h = slopes[reg-1]*(i - i_ref)*self.dx + self.y_peaks[reg-1][1]
                hs[i] = [h,h]
                
            for j in range(self.Ny):
                y = self.ys[j]
                
                if j == 0: #lower boundary
                    grid[j,i] = 0

                   
                elif i == 0: #inlet boundary 
                    if y <= self.y_peaks[0][1]: 
                        grid[j,i] = 0
                    else:
                        grid[j,i] = -1
                    
                elif i == self.Nx-1:# outlet boundary
                    if y <= self.y_peaks[-1][0]: 
                        grid[j,i] = 0
                    else:
                        grid[j,i] = -1
                
                elif i == i_ref: # pwl region change              
                    if math.isclose(y, h_left) or math.isclose(y, h_right): # true boundary point at region change
                        grid[j,i] = 0
                        
                    elif h_left < h_right: 
                        
                        if h_left < y and y < h_right: # x=h(y) vertical boundary
                            grid[j,i] = 0

                        elif y < h_left: # below vert boundary (interior)
                            grid[j,i] = 1
                            
                        else:             # below vert boundary(exterior)
                            grid[j,i] = -1
                
                    else:
                        if h_left > y and y > h_right: # x=h(y) vertical boundary
                            grid[j,i] = 0
                            
                        elif y < h_right: # below vert boundary (interior)
                            grid[j,i] = 1
                            
                        else:              # below vert boundary (exterior)
                            grid[j,i] = -1

                else:
                    if math.isclose(y,h): # true boundary point not at region change (from dx | slope)
                        grid[j,i] = 0

                    elif y < h:            # below boundary (interior)
                        grid[j,i] = 1
                    
                    else:                   # below boundary (exterior)
                        grid[j,i] = -1
        
        self.space = grid
        self.hs = hs

#------------------------------------------------------------------------------
# Boundary conditions on stream and velocity
#------------------------------------------------------------------------------
   
    def streamInlet(self, j):

        y = self.y0 + j * self.dy
        h = self.y0 + self.H_in
        
        if y <= h:
            u_term = self.U * y * (h-y)**2 / h**2
            dp_term = self.flux * y**2 * (3*h - 2*y) / h**3
            psi = u_term + dp_term 
            return psi
        else:
            return 0
    
    def velInlet(self, j):
        y = self.y0 + j * self.dy
        h = self.y0 + self.H_in
        
        if y <= h:
            u = self.dp_in * (y**2 - h*y)/2 + self.U * (h-y)/h
            return u
        else: 
            return 0
        
#------------------------------------------------------------------------------
# Boundary interpolation
#------------------------------------------------------------------------------
    def interp(self, scale, v_opp, v_bdry=0):
        v_nbr = v_bdry + (v_bdry - v_opp)*scale

        return v_nbr

    def scale_N(self, i,j): #N: (s=i, t=j+1)
        # x = self.xs[i]     #int
        # y = self.ys[j] 
        h = self.hs[i][0]
        
        # x_N = x           #ext
        y_N = self.ys[j+1] 
        
        # x_S = x           #opp
        y_S = self.ys[j-1] 
        
        # x_bdry = x
        y_bdry = h 
        
        # l1 = np.sqrt((x_N - x_bdry)**2 + (y_N - y_bdry)**2)
        l1 = np.abs(y_N-y_bdry)
        # l2 = np.sqrt((x_S - x_bdry)**2 + (y_S - y_bdry)**2)
        l2 = np.abs(y_S-y_bdry) 

        scale = l1/l2        

        return scale
    


    def scale_E(self, i,j): #E: (s=i+1, t=j)
        x = self.xs[i]
        y = self.ys[j]
        h_l = self.hs[i][0]
        h_r = self.hs[i][1]
        
        x_W = self.xs[i-1]
        # y_W = y
        
        x_E = self.xs[i+1]
        # y_E = y
        h_E_l = self.hs[i+1][0]
        
        G = (h_E_l - h_r)/(x_E - x)
        x_bdry =  x + (y - h_l)/G
        # y_bdry = y
        
        # l1 = np.sqrt((x_E - x_bdry)**2 + (y_E - y_bdry)**2)
        l1 = np.abs(x_E-x_bdry)
        # l2 = np.sqrt((x_W - x_bdry)**2 + (y_W - y_bdry)**2)
        l2 = np.abs(x_W-x_bdry) 

        scale = l1/l2     

        return scale
    
    def scale_W(self, i,j): #s=i-1, t=j
        x = self.xs[i]
        y = self.ys[j]
        h_l = self.hs[i][0]
        h_r = self.hs[i][1]
        
        x_E = self.xs[i+1]
        # y_W = y
        
        x_W = self.xs[i-1]
        # y_W = y
        h_W_r = self.hs[i-1][1]
        
        G = (h_W_r - h_l)/(x_W - x)
        x_bdry =  x + (y - h_r)/G

        # l1 = np.sqrt((x_W - x_bdry)**2 + (y_W - y_bdry)**2)
        l1 = np.abs(x_W-x_bdry)
        # l2 = np.sqrt((x_E - x_bdry)**2 + (y_E - y_bdry)**2)
        l2 = np.abs(x_E-x_bdry) 
        
        scale = l1/l2 

        return scale


    
    def scale_NE(self, i,j): #s=i+1, t=j+1
        x = self.xs[i]
        y = self.ys[j]
        h_r = self.hs[i][1]
        h_l = self.hs[i][0]
        
        x_E = self.xs[i+1]
        y_N = self.ys[j+1]
        h_E_l = self.hs[i+1][0]
        
        x_W = self.xs[i-1]
        y_S = self.ys[j-1]
        
        G = (h_E_l - h_r)/(x_E - x)
        x_bdry = x + (y - h_l)/(G-1)
        y_bdry = y + x_bdry - x

        l1 = np.sqrt((x_E-x_bdry)**2 + (y_N-y_bdry)**2)
        l2 = np.sqrt((x_W-x_bdry)**2 + (y_S-y_bdry)**2)

        scale = l1/l2 

        return scale
        
    def scale_SW(self, i,j): #SW: (i-1, j-1)
        x = self.xs[i]
        y = self.ys[j]
        h_r = self.hs[i][1]
        h_l = self.hs[i][0]
    
        x_W = self.xs[i-1]
        y_S = self.ys[j-1]
        h_W_r = self.hs[i-1][1]
        
        x_E = self.xs[i+1]
        y_N = self.ys[j+1]
        
        G = (h_W_r - h_l)/(x_W - x)
        x_bdry = x + (y - h_r)/(G-1)
        y_bdry = y + x_bdry - x

        l1 = np.sqrt((x_W-x_bdry)**2 + (y_S-y_bdry)**2)
        l2 = np.sqrt((x_E-x_bdry)**2 + (y_N-y_bdry)**2)

        scale = l1/l2 
   
        return scale
    
    def scale_NW(self, i,j): #NW: (s=i-1, t=j+1)
        x = self.xs[i]
        y = self.ys[j]
        h_r = self.hs[i][1]
        h_l = self.hs[i][0]
        
        x_W = self.xs[i-1]
        y_N = self.ys[j+1]
        h_W_r = self.hs[i-1][1]

        
        x_E = self.xs[i+1]
        y_S = self.ys[j-1]
        # h_SE = self.hs[i+1][0]
        
        G = (h_W_r - h_l)/(x_W - x)
        x_bdry = x + (y - h_r)/(G+1)
        y_bdry = y - x_bdry + x


        l1 = np.sqrt((x_W-x_bdry)**2 + (y_N-y_bdry)**2)
        l2 = np.sqrt((x_E-x_bdry)**2 + (y_S-y_bdry)**2)

        scale = l1/l2   
        # if scale > 1:
        #     print('nw', scale)
            
        return scale
    
    def scale_SE(self, i,j): #SE: (s=i+1, t=j-1)
        x = self.xs[i]
        y = self.ys[j]
        h_r = self.hs[i][1]
        h_l = self.hs[i][0]
        
        x_E = self.xs[i+1]
        y_S = self.ys[j-1]
        h_E_l = self.hs[i+1][0]
        
        x_W = self.xs[i-1]
        y_N = self.ys[j+1]

        G = (h_E_l - h_r)/(x_E - x)
        x_bdry = x + (y - h_l)/(G+1)
        y_bdry = y - x_bdry + x

        l1 = np.sqrt((x_E-x_bdry)**2 + (y_S-y_bdry)**2)
        l2 = np.sqrt((x_W-x_bdry)**2 + (y_N-y_bdry)**2)
        

        scale = l1/l2

        return scale
#------------------------------------------------------------------------------















