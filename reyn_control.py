# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:05:05 2024

@author: sarah
"""
import numpy as np
import graphics

import time

import reyn_boundary as bc


log_linthresh=1e-8  
        
class Reyn_Solution:
    def __init__(self, height, BC, pressure, velocity, solver_str, t=None):
        self.height = height
        self.BC = BC
        self.pressure = pressure
        self.velocity = velocity
        self.solver_str = solver_str
        
        if time is not None:
            self.time=t
        
        if isinstance(self.BC, bc.Mixed):
            self.dP = pressure.get_dP(height)
            self.Q = self.BC.Q
            
        elif isinstance(self.BC, bc.Fixed):
            self.dP = self.BC.dP
            self.Q = velocity.get_flux(height)
            
    def p_plot(self, scaled=False, zoom=False):

        if self.Q == 0: #cavity flow 
            paramstr = "$Q=%.2f$, $U=%.1f$"%(self.Q, self.BC.U)
        else:
            paramstr = "$Q=%.2f$, $U=%.1f$, $\Delta p=%.2f$"%(self.Q, self.BC.U, self.dP)
        
        p_title = self.solver_str +'\n' + paramstr
        p_labels = ["$p$", "$x$","$y$"]
       
        if zoom: 
            xs_zoom, ys_zoom = graphics.grid_zoom_1D(self.height.xs, self.height.ys, self.height )
            p_zoom = graphics.grid_zoom_2D(self.pressure.ps_2D, self.height )
            graphics.plot_contour_mesh(p_zoom, xs_zoom, ys_zoom, p_title, p_labels)
        else:
            graphics.plot_contour_mesh(self.pressure.ps_2D, self.height.xs, self.height.ys, p_title, p_labels)
        
    
    def v_plot(self, scaled=False, zoom=False,  inc=False, uv=False):
                
    
        if self.Q == 0: #cavity flow 
            paramstr = "$Q=%.2f$, $U=%.1f$"%(self.Q, self.BC.U)
        else:
            paramstr = "$Q=%.2f$, $U=%.1f$, $\Delta p=%.2f$"%(self.Q, self.BC.U, self.dP)
        
        v_title = self.solver_str+ '\n' + paramstr
        v_ax_labels =  ['$|(  u,  v)|_2$','$x$', '$y$']  
        uv_mag = np.sqrt((self.velocity.u)**2 + (self.velocity.v)**2)
        

        if zoom:
            xs_zoom, ys_zoom = graphics.grid_zoom_1D(self.height.xs, self.height.ys, self.height)
            u_2D_zoom = graphics.grid_zoom_2D(self.velocity.u, self.height)
            v_2D_zoom = graphics.grid_zoom_2D(self.velocity.v, self.height)
            uv_mag_zoom = graphics.grid_zoom_2D(uv_mag, self.height)
     
            graphics.plot_stream_heat(u_2D_zoom, v_2D_zoom, xs_zoom, ys_zoom, uv_mag_zoom, v_title, v_ax_labels) 
        else:
        
           graphics.plot_stream_heat(self.velocity.u, self.velocity.v, self.height.xs, self.height.ys, uv_mag, v_title, v_ax_labels)

        if uv:
            
            graphics.plot_contour_mesh(self.velocity.u, self.height.xs, self.height.ys, 'u', ['$u$', '$x$', '$y$'], -3, 3)
            graphics.plot_contour_mesh(self.velocity.v, self.height.xs, self.height.ys, 'v', ['$v$', '$x$', '$y$'], -3, 3)
            
        if inc:
            
            inc = self.velocity.make_inc(self.height)
            graphics.plot_contour_mesh(inc, self.height.xs, self.height.ys, 'incompressibility', ['$u_x+v_y$', '$x$', '$y$'], -1, 1)
            
            
            qs = self.velocity.get_flux_all(self.height)
           
            graphics.plot_2D(qs, self.height.xs, 'flux $\mathcal{Q} = q(x) =\int_0^{h(x)} u(x,y) dy$', ['$x$', '$q(x)=\mathcal{Q}$'])
