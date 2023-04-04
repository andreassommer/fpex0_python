"""
Self-written integrators for the solution of the semi-discretized Fokker-Planck equation.

31/03
Integrator with updating / moving grid.

Requirements.
- full control of solve_ivp
- take custom update functions for the grid
"""
# from scipy.integrate import solve_ivp
import numpy as np

"""
def vargrid_integrate(grid, update, fix_steps, integrator=solve_ivp, **kwargs):
    # DEV/NOTE: DRAFT
    # The solution must be stored every time the grid is changing.
    # A difficulty lies in keeping track of the changing grid.
    # In some way, a functionality is needed that stores all grids and
    # provides them when the residual is calculated.

    # Parameters could be read from FPEX0setup.
    # But: Function signature alterations would be necessary.
    # Also, big amounts of data might be passed every time.
    # Depends how Python passes class instances.
    t0tf = grid.gridTdot[[0,-1]] # extract
    first_step = None
    while tf_k <= t0tf[-1]:
        sol_k = integrator(first_step=first_step, max_step=fix_steps, **kwargs)
        
        # parameters for restart
        tf_k = sol_k.t[-1]
        first_step = sol_k.t[-1] - sol_k.t[-2]
        kwargs["t_span"] = [tf_k, t0tf[-1]]
        
        # do something
        update(grid)

        # interpolate
        # -- scipy.interpolate
        kwargs["y0"] = None # new initial value


    return sol_k
"""

def equidistribute_grid(Grid, U, alpha):
    """
    Determines on the basis of a given discrete solution U at time t an x-grid distributing
    the integral of the monitor function M(x,t) = (alpha + |u_xx(x,t)|)^1/2 over its intervals.
    This concept is known under the name of "equidistribution principle", 
    see e.g. Blom et al. on page 196, cited below.  
    
    ## Takes
    **Grid**  
    Grid object.

    **U**  
    Discrete solution calculate equidistributed grid on.

    **alpha**  
    Grid regularization parameter. See Blom et al., 1988.

    
    ## Returns
    **Grid**  
    x-grid (list) equidistributing the monitor function.

    ### Sources
    [1] Blom et al.: "On Simple Moving Grid Methods for One-Dimensional
        Evolutionary Partial Differential Equations"  
        in Journal of Computational Physics 74 (1988), 191-213.
    [2] Senz-Serna and Christie: "A Simple Adaptive Technique for Nonlinear Wave Problems"  
        in Journal of Computational Physics 67 (1986), 348-360.
    """
    # DEV/NOTE: might even customize monitor function

    # Algorithm from Senz-Serna and Christie (p. 351) with modifications by Blom et al.
    # It could be used any appropriate monitor function, for now we sticked to 
    # M(x,t) = (alpha + |u_xx(x,t)|)^1/2.
    
    # (1).
    J = len(Grid.gridT)-1
    assert(J==(len(U)-1))
    x = Grid.gridT
    h = Grid.h

    # (2).
    S_0 = 0
    S = [S_0] # store in S

    # (3). Calculate integral of the monitor by midpoint rule. Finite differences by Blom et al.
    # first interval
    FD = 2*( (U[2]-U[0])/(x[0]-x[2]) - (U[0]-U[1])/(x[0]-x[1]) )/(x[2]-x[0]) # u_xx((x_0+x_1)/2)
    r = np.sqrt( alpha + abs(FD) )
    S.append(S[-1]+r)
    # inner intervals
    for j in range(1, J-1):
        FD = ((U[j+2]-U[j])/(x[j+2]-x[j]) - (U[j+1]-U[j-1])/(x[j+1]-x[j-1]))/(x[j+1]-x[j]) 
        r = np.sqrt( alpha + abs(FD) )
        S.append(S[-1]+r)
    # last interval
    r = np.sqrt( alpha + abs(2*( (U[J]-U[J-1])/(x[J]-x[J-1]) - (U[J]-U[J-2])/(x[J]-x[J-2]) )/(x[J-1]-x[J-2])) )
    S.append(S[-1]+r)

    # (4). We skip assigning U^n+1.
    delta = S[-1]/J
    x_new = [x[0]]
    k=1

    # (5). 
    for j in range(1, J):
        B = j*delta
        # (6).
        while B <= S[k]:
            # (7).
            k+=1
        # (8). Calculate next grid point.
        # Note: The integration of the monitor function yields a polygon too, 
        #       so we can use the same formula.
        x_j_new = x[k-1] + (B-S[k-1])*h[k-1]/(S[k]-S[k-1]) # (we defined h differently)
        x_new.append(x_j_new)

        # (9). We don't assign U^n+1, thus skip this step.
    # skipped x_J in (4)
    x_new.append(x[-1])
    
    return x_new