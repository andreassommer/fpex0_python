import sys
import traceback

import numpy as np
from scipy import optimize
import time


def simulate(FPEX0setup, pvec, odeoptions={}, sensitivities=False):
    """
    Simulates Fokker-Planck with specified parameters for FP drift, diffusion, and initial function.

    ## Takes
    **FPEX0setup**
    <br> An FPEX0 setup configuration (`fpex0.setup.Setup`).

    **pvec**
    <br> Vector of parameters, coupled to FPEX0setup.
    <br> It must have the same scheme as FPEX0setup.Parameters.p_0, 
    then simulate() knows how to extract the parameters correctly.
    Usually this is ensured by the optimizer, who got the initial parameters and changes them through
    optimization steps.

    **odeoptions**
    <br> kwargs passed to the scipy solve_ivp ODE solver 
    (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html). 

    **method**
    <br> The method scipy solve_ivp should use.
    <br> Default is BDF, an an implicit multistep method.


    ## Returns
    **solution**
    <br> A scipy solve_ivp bunch object. Important parameters are
    - t: time points
    - y: corresponding values
    - sol: A (here) callable solution

    See link above for details.

    ### Comments
    The function will not let you set dense_output = False in odeoptions.
    """
    # extract parameters
    p_FPdrift     = FPEX0setup.Parameters.extract_p_FPdrift(pvec)
    p_FPdiffusion = FPEX0setup.Parameters.extract_p_FPdiffusion(pvec)
    p_IC          = FPEX0setup.Parameters.extract_p_iniDist(pvec)

    # evaluate initial distribution
    gridT         = FPEX0setup.Grid.gridT
    u0            = FPEX0setup.IniDist.f(gridT, p_IC)

    # retrieve "time" horizon for integrator
    t0tf          = FPEX0setup.Grid.gridTdot[[0, -1]]
    
    # generate right hand side, jacobian
    # insert control flow to compute either the nominal solution, the jacobian or both
    FPrhs         = FPEX0setup.make_rhsFcn(p_FPdrift, p_FPdiffusion, sensitivities)
    FPjac         = FPEX0setup.make_jacFcn(p_FPdrift, p_FPdiffusion, sensitivities, FPEX0setup.Integration.banded_jac)

    # setup integrator and update options, jacobian therein
    integrator  = FPEX0setup.Integration.integrator
    method      = FPEX0setup.Integration.method
    odeoptions  = FPEX0setup.Integration.updateOptions(odeoptions)
    odeoptions["dense_output"] = True   # dense_output for evaluable solution
    odeoptions  = FPEX0setup.Integration.updateJacobian(FPjac)
    
    # if sensivities are requested, the initial condition u0 has to be augmented
    if sensitivities:
        dy0dp = FPEX0setup.IniDist.dfdp(gridT, p_IC) # we received a gradient (matrix as we have requested several points)
        dy0dp = np.transpose(dy0dp)                  # gradients --> jacobian
        nu = len(gridT)     # size of the system's nominal part
        np_FP = len(p_FPdrift) + len(p_FPdiffusion)
        Gp0 = np.column_stack( (np.zeros((nu,np_FP)), dy0dp) ) # assemble jacobian. is zero for the params y0 does not depend on (FP params)
        u0 = np.append(u0, Gp0.flatten(order='F'))             # order='F': flatten in column-major order
    
    # start integration
    try:
        sTime = time.time()
        solution = integrator(FPrhs, t0tf, u0, method,**odeoptions)
        # DEV/NOTE: solve_ivp takes flag "vectorized". we should use it
        duration = time.time() - sTime
        print(f"Simulate: {duration:.3f}s")
    except:
        exception = sys.exc_info()
        traceback.print_exception(*exception)
        print('Integration failed!')
        duration = None
        solution = None
    
    return solution


def fit(FPEX0setup, optimizer='lsq'):
    """
    Fits the Fokker-Planck simulation to the given measurements as an optimization of the parameters
    for drift, diffusion and the initial distribution.

    ## Takes
    **FPEX0setup**
    <br> An FPEX0 setup configuration (`fpex0.setup.Setup`).

    **optimizer**
    <br> The optimizer that should be used. So far only least squares is implemented.


    ## Returns
    **result** 
    <br> A scipy.optimize.OptimizeResult object, with fit.x holding the parameter vector found.

    """
    # retrieve parameter values, bounds, indices
    p_0  = FPEX0setup.Parameters.p0
    p_lb = FPEX0setup.Parameters.p_lb
    p_ub = FPEX0setup.Parameters.p_ub
    
    # set function that computes the residual vector
    resvecfun = lambda p: residual(FPEX0setup, p)
    jac       = lambda p: jacobian(FPEX0setup, p)

    # optimization
    print(f'Running {optimizer}.\n')
    
    if optimizer.lower() == 'lsq':
        lsq_opts = {}
        # set options
        if FPEX0setup.Integration.vde:
            lsq_opts["jac"] = jac
        else:
            lsq_opts["jac"] = "3-point"

        lsq_opts["max_nfev"]                    = 100000    # max function evaluations
                                                            # tolerances
        lsq_opts["xtol"]                        = 1e-10      # x-step
        lsq_opts["ftol"]                        = 1e-20     # function-step
        lsq_opts["gtol"]                        = 1e-4       # norm of gradient, quite high, but okay for FD

        lsq_opts["x_scale"]                     = 'jac'     # let scipy set scale with jacobian
        lsq_opts["verbose"]                     = 2         # give detailed progress information
        result = optimize.least_squares(resvecfun, p_0, bounds=(p_lb, p_ub), **lsq_opts)
        
        return result

    else:
        raise ValueError("Your specified optimizer is not yet implemented. \nYou are welcomed to contact us by mail if interested in contribution!")


def residual(FPEX0setup, p_all):
    """
    Calculates the residual vector of measurements and simulation values, i.e. measVals - simVals
    at suitable points.
    
    
    ## Takes
    **FPEX0setup**
    <br> An FPEX0 setup configuration (`fpex0.setup.Setup`).

    **p_all**
    <br> Vector of parameters, coupled to FPEX0setup.
    <br> It must have the same scheme as FPEX0setup.Parameters.p_0, 
    then residual() knows how to use the parameters correctly.
    Usually this is ensured by the optimizer, who got the initial parameters and changes them through
    optimization steps.

    ## Returns
    **resvec**
    <br> Residual vector as described above.

    """
    measurements = FPEX0setup.Measurements
    meas_count   = len(measurements.rates)
    meas_values  = measurements.values
    meas_T       = measurements.temperatures
    meas_rates   = measurements.rates

    grid_T       = FPEX0setup.Grid.gridT

    # simulate
    # sTime = time.time()
    sol = simulate(FPEX0setup, p_all)
    # duration = time.time() - sTime
    # evaluate at measurement rates
    simdata = sol.sol(meas_rates)

    # residual
    resvec = np.array([])
    for k in range(meas_count):
        # select grid points matching to measurements
        _, idxGrid, idxMeas = np.intersect1d(grid_T, meas_T[k], assume_unique=True, return_indices=True)
        if len(idxGrid) != len(meas_T[k]):
            raise ValueError("Grid does not fit.")
        # get corresponding measurement and simulation data
        measVals    = meas_values[k][idxMeas]
        simVals     = simdata[idxGrid, k]
        resvec   = np.append(resvec, simVals - measVals)

    return resvec

def jacobian(FPEX0setup, p_all):
    """
    Docs.
    """
    # DEV/NOTE: We let the jacobian be computed by simulate. 
    #           By the use of flags we request the solution, it's jacobian or both.
    measurements = FPEX0setup.Measurements
    meas_count   = len(measurements.rates)
    meas_values  = measurements.values
    meas_T       = measurements.temperatures
    meas_rates   = measurements.rates

    grid_T       = FPEX0setup.Grid.gridT

    # simulate
    sol = simulate(FPEX0setup, p_all, sensitivities=True)
    # evaluate at measurement rates
    simdata = sol.sol(meas_rates)
    # extract sensitivities
    n_p = len(p_all)
    nu = len(sol.sol(grid_T[0]))//(n_p+1) # one nominal solution + for every parameter the sensivity (dy/dp)
    sensitivitydata = simdata[nu:,:]

    # residual
    jac = []
    for i, p in enumerate(p_all):
        dres_dp = np.array([])
        for k in range(meas_count):
            # select grid points matching to measurements
            _, idxGrid, _idxMeas = np.intersect1d(grid_T, meas_T[k], assume_unique=True, return_indices=True)
            if len(idxGrid) != len(meas_T[k]):
                raise ValueError("Grid does not fit.")
            
            # get sensivity data corresponding to measurements
            sensVals = sensitivitydata[idxGrid+i*nu, k]
            dres_dp = np.append(dres_dp, sensVals)

        jac.append(dres_dp)
    jac = np.column_stack(jac)
    
    return jac