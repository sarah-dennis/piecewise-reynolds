# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:49:45 2024

@author: sarah
"""
import numpy as np
import time
import readwrite as rw
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

def run_spLU(ex, u, v, old_psi, iters, past_iters, error_mod, write_mod, err_tol):
    print("file: %s"%ex.filestr)
    t0 = time.time() 
    M = Dpsi_cscmatrixBuild(ex)
    LU = splu(M)
    t_k0 = time.time()
    print("N=%d  build-time %.2fs"%(ex.N, t_k0-t0))
    max_err = 1
    for k in range(iters): 

        u, v = uv_approx(ex, u, v, old_psi)
        
        rhs = update_rhs(ex, u, v, old_psi)
    
        psi = LU.solve(rhs)


        if k % error_mod == 0: 
            max_err = np.max(np.abs(old_psi - psi))
            print("    k=%d  error: %.2e"%(k, max_err))
            if max_err < err_tol:
                break

        if k % write_mod == 0:
            rw.write_stokes(ex, u, v, psi, k+1+past_iters)
            
        old_psi = psi
    t_kf = time.time()
    print("N=%d  cnvg-error:%.2e   avg-iter-time:%.2fs"%(ex.N, max_err, (t_kf-t_k0)/(k+1)))
    return u, v, psi

class DPsi_Mat():
    def __init__(self, Ny, Nx):
        self.length = Ny * Nx * 9
        self.coefs = np.zeros(self.length)
        self.row = np.zeros(self.length)
        self.col = np.zeros(self.length)
        self.s = 0
        
    def append(self, k, nbr_k, coef):
        self.row[self.s] = k
        self.col[self.s] = nbr_k
        self.coefs[self.s] = coef
        self.s += 1
        
        
def Dpsi_cscmatrixBuild(ex):
    m = ex.Ny
    n = ex.Nx

    space = ex.space
    mat = DPsi_Mat(m, n)

    for k in range(m*n):
    
        i = k % n
        j = k // n
            
        # exterior -> identity row = 0
        # boundary -> identiy row = b.c.
        if space[j,i] != 1:
            mat.append(k, k, 1)
            
        # interior 
        # [... 1 -8  1 ... -8  28 -8 ... 1 -8  1 ...]
        else: #append(row:k, col:nbr(k), coef)
            mat.append(k, k, 28)
            
            mat.append(k, j*n + i-1, -8)
            mat.append(k, j*n + i+1, -8)
            mat.append(k, (j-1)*n + i, -8)
            mat.append(k, (j+1)*n + i, -8)
                
            mat.append(k, (j-1)*n + i-1, 1)  
            mat.append(k, (j-1)*n + i+1, 1)
            mat.append(k, (j+1)*n + i-1, 1)
            mat.append(k, (j+1)*n + i+1, 1)

    csc_mat = csc_matrix((mat.coefs, (mat.row, mat.col)), (m*n, m*n))

    return csc_mat


# -----------------------------------------------------------------------------
# rhs = c0*A + c1*(B - C) + Dpsi_bc

# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

def update_rhs(ex, u, v, psi): #
    
    n = ex.Nx 
    m = ex.Ny
    space = ex.space
    rhs = np.zeros(n*m)

    c0 = 3 * ex.dx
    c1 = 0.5 * ex.dx**2 * ex.Re
    
    for k in range(n*m):
        i = k % n
        j = k // n
            
        if i == 0 and space[j,i] == 0: # inlet: psi(y,x0) ~ Q
            rhs[k] = ex.streamInlet(j)
        
        elif i == n-1 and space[j,i] == 0: # outlet: dx{psi(y,x)} = 0
            rhs[k] = psi[j*n + i-1]   
            
        elif j == 0: # lower boundary y=0 
            rhs[k] = 0
            
        elif space[j,i] == 0: # upper boundary  y = h(x)
            rhs[k] = ex.flux
        
        elif space[j,i] == -1: #exterior
            rhs[k] = 0
        
        else: # interior
            # (u,v) at 9 point stencil
            u_k = u[k]
            v_k = v[k]
                   
            #possible exterior nbrs 
            dpsi_bc = 0
            k_N = (j+1)*n + i
            k_E = j*n + i+1
            k_W = j*n + i-1 
            k_S = (j-1)*n + i 
            k_NW = (j+1)*n + i-1
            k_NE = (j+1)*n + i+1
            k_SW = (j-1)*n + i-1
            k_SE = (j-1)*n + i+1
            
            # South (i, j-1)
            u_S = u[k_S]
            v_S = v[k_S]

            # North (i, j+1) 
            if ex.space[j+1,i] == -1: #N:
                scale_N = ex.scale_N(i,j)
                dpsi_bc += -8 * ex.interp(scale_N, psi[k_S], v_bdry= ex.flux)
                u_N = ex.interp(scale_N, u[k_S])
                v_N = ex.interp(scale_N, v[k_S])

            else:
                u_N = u[k_N]
                v_N = v[k_N]
            
            # East (i+1, j)                
            if ex.space[j,i+1] == -1: #E:

                scale_E = ex.scale_E(i,j)
                dpsi_bc += -8 * ex.interp(scale_E, psi[k_W], v_bdry= ex.flux)
                u_E = ex.interp(scale_E, u[k_W])
                v_E = ex.interp(scale_E, v[k_W])
                
            else:
                u_E = u[k_E]
                v_E = v[k_E]
        
            # West (i-1, j)          
            if ex.space[j,i-1] == -1: #W:
                scale_W = ex.scale_W(i,j)
                dpsi_bc += -8 * ex.interp(scale_W, psi[k_E], v_bdry= ex.flux)
                u_W = ex.interp(scale_W, u[k_E])
                v_W = ex.interp(scale_W, v[k_E])
            else:
                u_W = u[k_W]
                v_W = v[k_W]

            if ex.space[j+1,i+1] == -1 : #NE:
                scale_NE = ex.scale_NE(i,j)
                dpsi_bc += ex.interp(scale_NE, psi[k_SW], v_bdry= ex.flux)

            if ex.space[j+1,i-1] == -1: #NW:
                scale_NW = ex.scale_NW(i,j)
                dpsi_bc += ex.interp(scale_NW, psi[k_SE], v_bdry= ex.flux)
                                                     
            if ex.space[j-1,i+1] == -1: #SE:
                scale_SE = ex.scale_SE(i,j)
                dpsi_bc += ex.interp(scale_SE, psi[k_NW], v_bdry= ex.flux)

            if ex.space[j-1,i-1] == -1: #SW:
                scale_SW = ex.scale_SW(i,j)
                dpsi_bc += ex.interp(scale_SW, psi[k_NE], v_bdry= ex.flux)

            A = u_S - u_N + v_E - v_W
            
            B = v_k * (u_E + u_W + u_N + u_S)
            
            C = u_k * (v_E + v_W + v_N + v_S)
 
    
            rhs[k] = c0 * A + c1 * (B - C) - dpsi_bc
                   
    return rhs
 
# Velocity <-> Stream update

# u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
# v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)

def uv_approx(ex, u, v, psi):
  

    n = ex.Nx
    m = ex.Ny

    U = ex.U
    
    c2 = 3/(4*ex.dx)
    c3 = 1/4

    for k in range(n*m):
        i = k % n
        j = k // n

        # y=0 moving boundary
        if j == 0:
            u[k] = U
            v[k] = 0  

        elif i == 0: #inlet: u(x,y) ~ Q,U
            u[k] = ex.velInlet(j)
            v[k] = 0 
                
        elif i == n-1: #outlet: dx{u(dx)} = 0

            u[k] = u[j*n+i-1]
            v[k] = 0 
                    
        # y=h(x) boundary & exterior
        elif (ex.space[j,i] != 1):
            u[k] = 0
            v[k] = 0 
                
        else: #interior
            # (u,v, psi) at 4 point stencil
            k_E = j*n + i+1
            k_W = j*n + i-1
            k_S = (j-1)*n + i
            k_N = (j+1)*n + i
            
            
            # South (i, j-1)
            if j-1 == 0:  # y=0 boundary
                u_S = ex.U
                psi_S = 0
                
            else:  # exterior 
                u_S = u[k_S]
                psi_S = psi[k_S]
            
            # North (i, j+1)
            if ex.space[j+1,i] == 0: # h(x) boundary
                u_N = 0
                psi_N = ex.flux
                
            elif ex.space[j+1,i] == -1: # exterior 
                scale_N = ex.scale_N(i,j)
                u_N = ex.interp(scale_N, u[k_S]) 
                psi_N = ex.interp(scale_N, psi[k_S], v_bdry= ex.flux)
                
            else:
                u_N = u[k_N]
                psi_N = psi[k_N]
                
            # East (i+1, j)
            if i+1 == n-1 and ex.space[j,i+1] == 0: # outlet
                v_E = 0

                psi_E = psi[k] 
                
                
            elif ex.space[j,i+1] == 0: # y=h(x) boundary
                v_E = 0
                psi_E = ex.flux 
                
            elif ex.space[j,i+1] == -1:  # exterior 

                scale_E = ex.scale_E(i,j)
                v_E = ex.interp(scale_E, v[k_W])
                psi_E = ex.interp(scale_E, psi[k_W], v_bdry= ex.flux)
            else:
                v_E = v[k_E]
                psi_E = psi[k_E] 

            # West (i-1, j)  
            if i-1 == 0 and ex.space[j,i-1] == 0:  # inlet
                v_W = 0
                psi_W = ex.streamInlet(j)
                
            elif ex.space[j,i-1] == 0: # y=h(x) boundary
                v_W = 0
                psi_W = ex.flux 
                
            elif ex.space[j,i-1] == -1:  # exterior 
                scale_W = ex.scale_W(i,j)
                v_W = ex.interp(scale_W, v[k_E])
                psi_W = ex.interp(scale_W, psi[k_E], v_bdry= ex.flux)
            else:
                v_W = v[k_W]
                psi_W = psi[k_W]
   

            
            u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)

    return u, v


        
        
        
        
        
        
        
        
        