# PIECEWISE REYNOLDS
Code for, 'FAST SOLVER FOR THE REYNOLDS EQUATION ON PIECEWISE LINEAR GEOMETRIES' (Sarah Dennis & Thomas Fai, 2026)

This project includes finite difference solvers for the Reynolds equation, the Stokes equation,
and piecewise exact solvers for the Reynolds equation.


#------------------------------------------------------------------------------
# REYNOLDS SOLVERS: FD, PWL, PWC 
#------------------------------------------------------------------------------

To solve the Reynolds equation use reyn_run.py

1. Choose an Example and set the geometry parameters (args)
    
    A.) backward facing step

          #--------------------------------------------------------------------
            Example = examples.BFS
            h_inlet = 2
            h_outlet = 1
            l_inlet = 8
            L_total = 16
            args =  [h_inlet, h_outlet, l_inlet, L_total]
          #--------------------------------------------------------------------
    
    B.) piecewise linear height (regularized BFS)
          
          #--------------------------------------------------------------------
            Example = examples.BFS_pwl
            H = 2
            h = 1
            delta = 1
            L = 16
            args = [H,h,L,delta]
          #--------------------------------------------------------------------
          
    C.) logistic step
    
          #--------------------------------------------------------------------
            Example = examples.Logistic
            delta = 8     
            h_in = 1       
            h_out = 2    
            L_total = 16     
            args = [ h_in, h_out, L_total, delta]
          #--------------------------------------------------------------------
          
    D.) sinusoidal slider
    
          #--------------------------------------------------------------------
            Example = examples.Sinusoid
            H = 1
            delta = 1/2
            k = 1 
            L = 2
            args = [H, delta, k, L]
          #--------------------------------------------------------------------

    --> See reyn_examples.py, reyn_heights.py, domain.py
    

2. Set the boundary conditions 

    - U: lower surface velocity
    - Q: flux 
    
          #--------------------------------------------------------------------    
            U = 0
            Q = 1
            BC = bc.Mixed(U,Q)
          #--------------------------------------------------------------------   
    
    --> See boundary.py

3. Initialize the Solver and grid scale (N)

    - Solver: class with Reynolds, ELT, PLT solve methods
    - N: grid scale (dx = 1/N)
    
          #--------------------------------------------------------------------    
            Solver = control.Reynolds_Solver(Example, BC, args)
            N = 160
          #--------------------------------------------------------------------    
          
    --> see reyn_control.py

4. Run solver
   
    A.) finite difference solver
    
          #--------------------------------------------------------------------    
            solution = solver.fd_solve(N)
          #--------------------------------------------------------------------
    
    B.) piecewise constant solver (PWC))
    
         #--------------------------------------------------------------------    
            solution = solver.pwc_schur_solve(N)
         #--------------------------------------------------------------------

    C.) piecewise linear solver (PWL)
   
        #--------------------------------------------------------------------    
           solution = solver.pwl_schur_solve(N)
        #--------------------------------------------------------------------
 
    --> see reyn_control.py, graphics.py, reyn_pressure.py, reyn_pressure_finDiff.py, 
        reyn_pressure_pwc.py, reyn_pressure_pwl.py, reyn_velocity.py
        

5. Plot the solution

    - set plots_on=True at the top of reyn_run

        #--------------------------------------------------------------------    
            if plots_on:
                solution.p_plot()
                solution.v_plot()
        #--------------------------------------------------------------------        

#------------------------------------------------------------------------------
# STOKES SOLVER
#------------------------------------------------------------------------------

To use the Stokes solver, use file stokes_run.py

1. Choose an Example and set the geometry parameters (args)
    
    - Examples are set identically to step 1. for reyn_run.py above
    
    --> See stokes_examples.py, stokes_heights.py, domain.py
    

2. Set the boundary conditions 

    - Boundary conditions are set identically to step 2. for reyn_run.py above
    
    --> See boundary.py

3. Initialize the Solver and grid scale (N)

    - Solver: class with iterative solve methods
    - N: grid scale (dx = 1/N)
    - Re: Reynolds number for stream-velocity Navier-Stokes equations (use Re=0 for Stokes solver)
    
          #--------------------------------------------------------------------    
            Solver = control.Stokes_Solver(Example, args, BC, Re=0) 
            
            N = 160
          #--------------------------------------------------------------------    
          
    --> see stokes_control.py

4. Run the iterative Stokes solver

    - The iterative Stokes solver computes the stream function and velocity; 
      in general, the pressure is not computed (unless using load(N) or load_plot(N))
    - The iterative solver will run until the error in the stream function between
      iterations falls below 10^-8, or until the maximum number of iterations has been reached
    - The solver checks this error and saves the current solution every 500 iterations
    - The parameters {err_tol, max_iters, error_mod, write_mod} can be updated in stokes_control.py
    - Solutions are written to ./stokes_examples/example_name/example_name_N

    A.) run a new solution at grid size N
    
          #--------------------------------------------------------------------    
            Solver.new_run(N)
          #--------------------------------------------------------------------    
    
    B.) load and run an existing solution at grid size N
    
          #--------------------------------------------------------------------    
            Solver.load_run(N)
          #--------------------------------------------------------------------      
          
    C.) run a new solution starting at grid size N,
        then interpolate the solution to grid size 2*N and run,
        repeat, ending with grid size (2^k)*N
    
          #--------------------------------------------------------------------
            k = 4    
            Solver.new_run_many(N, 2, k)
          #--------------------------------------------------------------------    
    
    D.) load and run an existing solution starting at grid size N,
        then interpolate the solution to grid size 2*N and run,
        repeat, ending with grid size (2^k)*N   

          #--------------------------------------------------------------------
            k = 4    
            Solver.load_run_new_many(N, 2, k)
          #--------------------------------------------------------------------    
        
    E.) load and run an existing solution starting at grid size N,
        then load and run an existing solution at grid size 2*N,
        repeat, ending with grid size (2^k)*N   

          #--------------------------------------------------------------------
            k = 4    
            Solver.load_run_many(N, 2, k)
          #--------------------------------------------------------------------    
    
    F.) load an existing solution at grid size N
        
    - this is the only method that computes the pressure and returns the solution

          #--------------------------------------------------------------------   
            p, u, v, stream = Solver.load(N)
          #--------------------------------------------------------------------    

    --> see stokes_control.py, stokes_solver.py, stokes_pressure.py, read_write.py

5. Plot the solution
    
    - the solution at grid size N must already exist (see step 4.)
    - the pressure is computed from the saved stream-velocity solution
    - the bool-flag zoom-on is set at the top of stokes_run.py
    
          #--------------------------------------------------------------------
            Solver.load_plot(N, zoom=zoom_on)
          #--------------------------------------------------------------------    
       
    --> see stokes_control.py,  graphics.py

 
#------------------------------------------------------------------------------
# TIMING & CONVERGENCE
#------------------------------------------------------------------------------

To test the timing and convergence of the Reynolds FD, PWC & PWL, methods
run the script reyn_piecewise_convg_timing.py

To see linear time performance for PWL, the 2D PRESSURE MUST BE TURNED OFF
    1. Go to reyn_pressure.py
    2. In the Reyn_Pressure class, set self.2D_ps = None
    Otherwise, self.2D_ps = self.make_2D_ps(height, ps_1D) sets p(x,y) = p(x) for 2D pressure plotting