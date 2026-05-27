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
zoom_on =not True         # plot a zoomed-in window, set location in graphics.py

#------------------------------------------------------------------------------
## Piecewise-linear examples 
#------------------------------------------------------------------------------

# Example = examples.BFS
# h_in=1
# h_out=2
# l_in=8
# l_out=8
# args =  [h_in, h_out, l_in, l_out]


# Example = examples.BFS_pwl
# H_in = 1
# H_out =2
# delta =2
# L =16
# args = [H_in,H_out,L,delta]


#------------------------------------------------------------------------------
## Smooth examples  
# ------------------------------------------------------------------------------
Example = examples.Sinusoid
H=1
delta = 1/2
k = 1 #* 2pi
L=2
args = [H, delta, k, L]

# Example = examples.Sinusoid_2
# H=1.5 # equilibrium
# delta = 1/3 # k even: hmax/H -1 = delta, k odd: 1 -hmax/H =delta
# k = 2 # period k * pi on length 2l
# l = 1 # half length of texture
# L=8 # half total length
# args = [H, delta, k, l, L]


# Example = examples.Logistic
# delta = 32 # max slope: delta*(H-h)/4
# H_in = 1      # inlet height
# H_out = 2      # outlet height
# L = 16     # total length
# args = [ H_in, H_out, L, delta]


#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

## U: velocity BC {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U = 0

#fixed pressure BC {p(x0,y)=-dP, p(xL,y)=0} 
dP = 50
BC = bc.Fixed(U,dP)

# mixed pressure BC {dp/dx (x0,y) ~ Q, p(xL,y)=0}
# Q = 1
# BC = bc.Mixed(U, Q)

#------------------------------------------------------------------------------

solver = solvers.Reynolds_Solver(Example, BC, args)

#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )

N =160
solution = solver.fd_lu_solve(N)

# solution = solver.fd_bicg_solve(N)

# solution = solver.pwc_solve(N)

# if __name__ == '__main__':
#     solution = solver.pwc_parallel_solve(N)

# solution = solver.pwl_solve(N)

# solution = solver.pwl_gmres_solve(N)

print('solve time: %.2fs'%solution.time)
if plots_on:
    solution.p_plot(zoom=zoom_on)
    solution.v_plot(zoom=zoom_on, uv=uv_on, inc=inc_on)
#------------------------------------------------------------------------------
