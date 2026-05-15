# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 11:23:51 2026

@author: sarah
"""

import reyn_solvers as solvers
import reyn_boundary as bc
import reyn_examples as examples
import graphics
import numpy as np

# TODO: 
# for linear timing with PWL method.... TURN OFF 2D PRESSURE
#   in          reyn_pressure.py, class Reyn_Pressure,
#   change      self.ps_2D = self.make_2D_ps(height, ps_1D)
#   to          self.ps_2D = None

#------------------------------------------------------------------------------
#TODO: set example

# Example = examples.BFS
# H=2
# h=1
# L=16
# args =  [H, h, L/2, L/2]


# Example = examples.BFS_pwl
# H_in = 2
# H_out =1
# delta =1/8
# L =16
# args = [H_in,H_out,L,delta]


Example = examples.Sinusoid
H=1
delta = 1/2
k = 1 #* 2pi
L=16
args = [H, delta, k, L]


#------------------------------------------------------------------------------
#TODO: Set boundary conditions

U = 2
#U = 0

Q = (U*H/2) * (1-(delta**2))/(1+(delta**2)/2) #sinuosoid Q for DP=0
#Q = 1


BC = bc.Mixed(U, Q)

solver = solvers.Reynolds_Solver(Example, BC, args)

#------------------------------------------------------------------------------

tests = 8  # N=2**k for k in k0+[0,tests)  
k_0=2           
                                                                                                                                                                                                          
Ns = np.zeros(tests)


dP_err_fd= np.zeros(tests)
dP_err_pwc= np.zeros(tests)
dP_err_pwl= np.zeros(tests)

l2_err_fd= np.zeros(tests)
l2_err_pwc= np.zeros(tests)
l2_err_pwl= np.zeros(tests)

pwc_schur_times= np.zeros(tests)
pwl_schur_times=np.zeros(tests)
fd_times= np.zeros(tests)


for k in range(tests):

    N = 2**(k_0+k)
    Ns[k]=N
    print(f'k={k+1:d} of {tests:d}, N={N:d}')
    
    fd_solution = solver.fd_solve(N)
    pwc_solution = solver.pwc_solve(N)
    pwl_solution = solver.pwl_solve(N)

    
#TODO: set exact solution
    exact_solution = solver.exact_sinusoid_sol(N)
    # exact_solution = pwl_solution


    dP_err_fd[k] = abs(fd_solution.dP-exact_solution.dP) 
    dP_err_pwl[k] = abs(pwl_solution.dP-exact_solution.dP)  
    dP_err_pwc[k] = abs(pwc_solution.dP-exact_solution.dP)
    
    l2_err_fd[k] = (sum((fd_solution.pressure.ps_1D-exact_solution.pressure.ps_1D)**2)/N)**(1/2)  
    l2_err_pwc[k] = (sum((pwc_solution.pressure.ps_1D-exact_solution.pressure.ps_1D)**2)/N)**(1/2)  
    l2_err_pwl[k] = (sum((pwl_solution.pressure.ps_1D-exact_solution.pressure.ps_1D)**2)/N)**(1/2)  

    fd_times[k] = fd_solution.time
    pwc_schur_times[k] = pwc_solution.time
    pwl_schur_times[k] = pwl_solution.time
    

#------------------------------------------------------------------------------        
# plot ||dp-dp*|| error
graphics.plot_log_multi([dP_err_fd, dP_err_pwc, dP_err_pwl], Ns, 'Convergence in $\Delta P$: Absolute error', ['FD', 'PWC', 'PWL'], ['$1/\Delta x$', 'error'], loc='upper', bigO_on=True)
# plot ||p-p*|| error
graphics.plot_log_multi([l2_err_fd, l2_err_pwc, l2_err_pwl], Ns, 'Convergence in $p(x)$: $l_2$ error', ['FD', 'PWC', 'PWL'], ['$1/\Delta x$', 'error'], loc='upper', bigO_on=True)

# plot FD, PWC, PWL times (linear vs quadratic)
graphics.plot_2D_multi([fd_times, pwc_schur_times, pwl_schur_times], Ns, 'Run Time', ['FD LU', 'PWC', 'PWL'], ['$1/\Delta x$', 'run time (s)'], loc='left')
# plot PWL times only (linear)
graphics.plot_2D(pwl_schur_times, Ns, 'Run Time for PWL', ['$1/\Delta x$', 'run time (s)'], color='forestgreen', marker='s')
