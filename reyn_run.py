# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""


import reyn_boundary as bc
import reyn_examples as examples
import reyn_solvers as solvers

#-------------------plotting---------------------------------------------------

plots_on = True
uv_on = False          # plot u(x,y) & v(x,y) & |(u,v)|
inc_on = False         # plot ux + vy =? 0
zoom_on =False        # plot a zoomed-in window, set location in reyn_solution.py

#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------

Example = examples.BFS
h_in=1 
h_out=2
l_in=8
l_out=8
args =  [h_in, h_out, l_in, l_out]


# Example = examples.BFS_pwl
# H_in = 2
# H_out =1
# delta = 1
# L =16
# args = [H_in,H_out,L,delta]


#------------------------------------------------------------------------------
## Smooth examples  
##      (finite difference solution only)
# ------------------------------------------------------------------------------
# Example = examples.Sinusoid
# H=1
# delta = 1/2
# k = 1 #* 2pi
# L=2
# args = [H, delta, k, L]

# Example = examples.Sinusoid_2
# H=1
# delta = 1/4
# k = 2 # period k * pi on length 2l
# l = 1 # half length of texture
# L=3 # half length total length
# args = [H, delta, k, l, L]


Example = examples.Logistic
delta = 8 # max slope: delta*(H-h)/4
H = 2     # outlet height
h = 1       # inlet height
L = 8     # total length
args = [ H, h, L, delta]


#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

## U: velocity BC {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U = 0

#fixed pressure BC {p(x0,y)=-dP, p(xL,y)=0} 
# dP = 8
# BC = bc.Fixed(U,dP)

# mixed pressure BC {dp/dx (x0,y) ~ Q, p(xL,y)=0}
Q = 1
BC = bc.Mixed(U, Q)

#------------------------------------------------------------------------------

solver = solvers.Reynolds_Solver(Example, BC, args)

#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )

N = 100
# solution = solver.fd_solve(N)
# 
# solution = solver.pwc_solve(N)

# if __name__ == '__main__':
#     solution = solver.pwc_parallel_solve(N)

solution = solver.pwl_solve(N)

# solution = solver.pwl_gmres_solve(N)

if plots_on:
    solution.p_plot(zoom=zoom_on)
    solution.v_plot(zoom=zoom_on, uv=uv_on, inc=inc_on)
#------------------------------------------------------------------------------
