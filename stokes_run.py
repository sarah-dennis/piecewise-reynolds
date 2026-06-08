# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples


zoom_on=  True    #set zoom location in graphics.py

U=0
Q=1

# U=1
# Q=0

Re=0


#------------------------------------------------------------------------------
h_in = 1
h_out = 2
l_in = 8
l_out=8
args = [h_in, h_out, l_in, l_out]
Example = examples.BFS

#------------------------------------------------------------------------------
# h_in=1
# h_out=2
# L=16    #l_in=l_out=(L-l_wedge)/2
# l_wedge= 2#1/8   
# args = [h_in, h_out, L, l_wedge]
# Example = examples.BFS_pwl

#------------------------------------------------------------------------------
# H_in=1
# H_out=2
# L=16
# delta = 32 #slope: delta*(H_in-H_out)/4
# args = [H_in, H_out, L, delta]
# Example = examples.Logistic

#------------------------------------------------------------------------------


# H=1.5       # equilibrium H
# delta = 1/3 # amplitude delta*H
# l = 1       # texture length 2l
# k = 2       # period k * pi on length 2l
# L=8         # total length 2L
# args = [H, delta, k, l, L]

# Example = examples.Sinusoid

#------------------------------------------------------------------------------
# l=7
# h=1
# H = 2
# args =  [h, H, h, l, 1.25, 0.75, l]
# Example = examples.TriSlider



#------------------------------------------------------------------------------
solver = control.Stokes_Solver(Example, args, U, Q, Re, max_iters=500000)                
  
N=160


# solver.new_run(N) 

# solver.load_run(N)

# solver.load_scale(N,2*N) 

# solver.load_copy(N, new_Example, new_args)

# solver.load_run_many(N, 2, 3)

# solver.new_run_many(N, 2, 4)  

# solver.load_run_new_many(N, 2,3)
   
# solver.load_plot(N, zoom=zoom_on)

# ------------------------------------------------------------------------------
# solver.compare(args, U, Q, Re, 10,[20,40,80],160)#],320)
solver.compare(args, U, Q, Re, 20,[40,80],160)#],320)
# 






