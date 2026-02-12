#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:59:46 2024

@author: sarahdennis
"""
import numpy as np
from scipy.interpolate import interpn

import graphics
import readwrite as rw
import convergence as cnvg
import stokes_pressure as pressure

#------------------------------------------------------------------------------
from stokes_solver import run_spLU
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
class Stokes_Solver:
    def __init__(self, Example, args, U, Q, Re, max_iters=50000):
        # domain 
        self.Example = Example
        self.args = args
        self.U = U
        self.Q = Q
        self.Re = Re
        
        # iterative solution args
        self.max_iters = max_iters
        self.write_mod = 500
        self.error_mod = 500
        self.err_tol = 1e-8

        
#------------------------------------------------------------------------------
    def new_run(self, N):

        ex = self.Example(self.args, self.U, self.Q, self.Re, N)
        
        u_init = np.zeros(ex.Nx * ex.Ny)
        v_init = np.zeros(ex.Nx * ex.Ny)
        psi_init = np.zeros(ex.Nx * ex.Ny)
        past_iters = 0
    
        u, v, psi = run_spLU(ex, u_init, v_init, psi_init, self.max_iters, past_iters, self.error_mod, self.write_mod, self.err_tol)
    
        rw.write_stokes(ex, u, v, psi, self.max_iters)
    
    def new_run_many(self, N_0, dN, many):
                
        
        self.new_run(N_0)
        N_load = N_0
        for k in range (1, many): 
            N = N_load*dN
            self.load_scale(N_load, N)
            self.load_run(N)
            N_load = N

    def load_run_new_many(self, N_0, dN, many):
        N_load = N_0
        self.load_run(N_load)
        for k in range (many): 
            N = N_load*dN
            self.load_scale(N_load, N)
            self.load_run(N)
            N_load = N

    def load_run_many(self, N_0, dN, many):
        N = N_0
        for k in range (many): 
            self.load_run(N)
            N *= dN 
                                                                                                                                                                                                                                                                       
    def load_run(self, N):                                
        ex = self.Example(self.args, self.U, self.Q, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx*ex.Ny)
        u, v, psi = run_spLU(ex, u, v, psi, self.max_iters, past_iters, self.error_mod, self.write_mod, self.err_tol)
        rw.write_stokes(ex, u, v, psi, self.max_iters+past_iters)
   
    def load_scale(self, N_load, N_scale):
        
        ex_load = self.Example(self.args, self.U, self.Q, self.Re, N_load)
        ex_scale = self.Example(self.args, self.U, self.Q, self.Re, N_scale)
        
        points_load = (ex_load.ys, ex_load.xs)
        
        u_load, v_load, psi_load, past_iters = rw.read_stokes(ex_load.filestr+".csv", ex_load.Ny*ex_load.Nx)
        u_load_2D = u_load.reshape((ex_load.Ny,ex_load.Nx))
        v_load_2D = v_load.reshape((ex_load.Ny,ex_load.Nx))
        psi_load_2D = psi_load.reshape((ex_load.Ny,ex_load.Nx))
    
        points_scale = np.meshgrid(ex_scale.ys, ex_scale.xs)
        
        u_scaled_2D = interpn(points_load, u_load_2D, tuple(points_scale))
        v_scaled_2D = interpn(points_load, v_load_2D, tuple(points_scale))
        psi_scaled_2D = interpn(points_load, psi_load_2D, tuple(points_scale))
    
        u_scaled = u_scaled_2D.ravel(order='F')
        v_scaled = v_scaled_2D.ravel(order='F')
        psi_scaled = psi_scaled_2D.ravel(order='F')
    
        rw.write_stokes(ex_scale, u_scaled, v_scaled, psi_scaled, 0)
    
    #intended for BFS -> wedge BFS
    def load_copy(self, N, new_Example, new_args):
        
        ex_load = self.Example(self.args, self.U, self.Q, self.Re, N)
        ex_new = new_Example(new_args, self.U, self.Q, self.Re, N)
        
        u_load, v_load, psi_load, past_iters = rw.read_stokes(ex_load.filestr+".csv", ex_load.Ny*ex_load.Nx)
    
        rw.write_stokes(ex_new, u_load, v_load, psi_load, 0)
    
    def load(self,N):
        ex = self.Example(self.args, self.U, self.Q, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx*ex.Ny)
        p = pressure.pressure(ex, u, v)
        dp = pressure.get_dp(ex, p)
        p = p.reshape((ex.Ny,ex.Nx))
        u = u.reshape((ex.Ny,ex.Nx))
        v = v.reshape((ex.Ny,ex.Nx))
        
        return p, u, v, dp
        
#------------------------------------------------------------------------------    
    def get_dP(self,N):
        ex = self.Example(self.args, self.U, self.Q, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)
        p = pressure.pressure(ex, u, v)
        dp = pressure.get_dp(ex, p) 
        return dp
    
    
    
#------------------------------------------------------------------------------
# Error
#------------------------------------------------------------------------------
    def compare(self, args, U, Q, Re, N_min, Ns, N_max,p_err=False): # grid convergence (multiple grid sizes of same example)
        
        l1_errs, l2_errs, inf_errs, cnvg_rates, ex_min = cnvg.stokes_cnvg_self(self.Example, args, U, Q, Re, N_min, Ns, N_max,p_err)
        title = "Iterative Grid Error in Stream $\psi$ at $N_{max}=%d$"%(N_max)
        ax_labels = ["$N$", "$||\psi _{N^{*}} - \psi_{N}||_p$"]
        leg_labels = ['$L^1$', '$L^2$','$L^\infty$']
        
        graphics.plot_log_multi([l1_errs, l2_errs, inf_errs], [N_min]+Ns, title, leg_labels, ax_labels,bigO_on=True,loc='lower' )

#------------------------------------------------------------------------------
# PLOTTING 
#------------------------------------------------------------------------------
    def load_plot(self, N,zoom=False):
        ex = self.Example(self.args, self.U, self.Q, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)

    
    # Grid domain
        xs = ex.xs
        ys = ex.ys
        
        
    # Zoom domain
        if zoom:
            xs_zoom, ys_zoom = graphics.grid_zoom_1D(xs, ys, ex)

    # Space plot: (debugging)
        # ax_labels = ['space', '$x$', '$y$']
        # if zoom:
        #     space_zoom = graphics.grid_zoom_2D(ex.space, ex)
        #     graphics.plot_contour_mesh(space_zoom, xs_zoom, ys_zoom,'space', ax_labels, vmin=-1, vmax=1)
        # else:
        #     graphics.plot_contour_mesh(ex.space, xs, ys, 'space', ax_labels, vmin=-1, vmax=1)
            
    # Pressure plot: 
        p = pressure.pressure(ex, u, v)
        dp = pressure.get_dp(ex, p) 
        
        p_2D = p.reshape((ex.Ny,ex.Nx))
        
        if self.Q == 0:
            dp_str = '' #cavity flow dp undefined
        else:
            dp_str = ', $\Delta p =%.2f$'%(dp)

    
        ax_labels_p = ['$p$', '$x$', '$y$']
        title_p = 'Stokes\n' + ex.spacestr + dp_str
    
        p_ma = np.ma.masked_where(ex.space!=1, p_2D)
       
        if zoom:
            p_zoom = graphics.grid_zoom_2D(p_ma, ex)  
            graphics.plot_contour_mesh(p_zoom, xs_zoom, ys_zoom, title_p, ax_labels_p)
        else:
            graphics.plot_contour_mesh(p_ma, xs, ys, title_p, ax_labels_p)
    
    #  Velocity plot: 
    
        ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
        
        title = 'Stokes\n' + ex.spacestr + dp_str
        ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
        
        u_2D = u.reshape((ex.Ny,ex.Nx))

        v_2D = v.reshape((ex.Ny,ex.Nx)) 

        uv_mag = np.sqrt(u_2D**2 + v_2D**2)
        

        if zoom:
            u_2D_zoom = graphics.grid_zoom_2D(u_2D, ex)
            v_2D_zoom = graphics.grid_zoom_2D(v_2D, ex)
            uv_mag_zoom = graphics.grid_zoom_2D(uv_mag, ex)
            graphics.plot_stream_heat(u_2D_zoom, v_2D_zoom, xs_zoom, ys_zoom, uv_mag_zoom, title, ax_labels)
        else:
            
            graphics.plot_stream_heat(u_2D, v_2D, xs, ys, uv_mag, title, ax_labels) 
        
    # Stream plot:
    
        # ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = -\psi_x$', '$x$', '$y$']
        # title = 'Stream $\psi(x,y)$ \n' + ex.spacestr + dp_str
        # stream_2D = psi.reshape((ex.Ny,ex.Nx))
        # stream_2D_ma = np.ma.masked_where(ex.space==-1, stream_2D)
        # if zoom:
        #     stream_2D_zoom = graphics.grid_zoom_2D(stream_2D_ma, ex)
        #     vmin = np.min(stream_2D_zoom)
        #     graphics.plot_contour_mesh(stream_2D_zoom, xs_zoom, ys_zoom, title, ax_labels, log_cmap=True,vmin=vmin, vmax=ex.flux)
        # else:
        #     graphics.plot_contour_mesh(stream_2D_ma, xs, ys, title, ax_labels, vmin=0, vmax=ex.flux)
        

