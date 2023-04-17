"""
Self-written integrators for the solution of the semi-discretized Fokker-Planck equation.
"""

import fpex0.setup as setup
from scipy.integrate import quad
from scipy.integrate import solve_ivp
import numpy as np


def equigrid_t0(FPEX0setup, m):
    # IN PROGRESS
    M0    = FPEX0setup.Integration.monitor
    x0xf  = FPEX0setup.Grid.gridT[[0,-1]]
    x0,xf = x0xf[0], x0xf[1]
    p0_IC = FPEX0setup.Parameters.p0_iniDist
    u_x   = lambda x: FPEX0setup.IniDist.dfdx( np.array([x]), np.array(p0_IC))
    u_xx  = lambda x: FPEX0setup.IniDist.dfdxx(np.array([x]), np.array(p0_IC))
    eta0 = quad(lambda x: M0(x, u_x, u_xx), x0, xf)[0]

    # integrate
    RHS = lambda t, x: eta0/M0(x, u_x, u_xx)
    sol = solve_ivp(fun=RHS, t_span=[0,1], y0=[x0xf[0]], dense_output=True, rtol=1e-5, atol=1e-8)
    x = sol.sol
    # construct grid
    t = np.linspace(0, 1, m)
    xgrid = x(t)[0]

    return xgrid


def equigrid(Grid, U, alpha):
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
    # DEV/TODO: make monitor function custom

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


def equigrid_odeint(update, FPEX0setup, gridpoints=1001, **kwargs):
    # DEV/NOTE: DRAFT
    # The solution must be stored every time the grid is changing.
    # A difficulty lies in keeping track of the changing grid.
    # In some way, a functionality is needed that stores all grids and
    # provides them when the residual is calculated.

    t0tf = FPEX0setup.Grid.gridTdot[[0,-1]] # extract
    integrator = FPEX0setup.Integration.integrator
    first_step = None
    fixed_steps = FPEX0setup.Integration.fixed_steps

    # an inital grid needs to be calculated first
    initial_gridT = equigrid_t0(FPEX0setup, gridpoints)
    initial_gridTdot = FPEX0setup.Grid.gridTdot
    initial_grid = setup.Grid(initial_gridT, initial_gridTdot, uniform=False)


    # other steps on updating grid
    while tf_k <= t0tf[-1]:
        sol_k = integrator(first_step=first_step, max_step=fixed_steps, **kwargs)
        
        # parameters for restart
        tf_k = sol_k.t[-1]
        first_step = sol_k.t[-1] - sol_k.t[-2]
        kwargs["t_span"] = [tf_k, t0tf[-1]]
        
        # do something
        grid = None # !DUMMY!
        update(grid)

        # interpolate
        # -- scipy.interpolate
        kwargs["y0"] = None # new initial value


    return sol_k

