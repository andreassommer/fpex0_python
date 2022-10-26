import importlib.resources as resources
import json
import numpy as np
import sympy

from fpex0.setup import Setup, Parameters, Grid, Integration
from fpex0 import FokkerPlanck
# from fpex0.default import defaultDiffusionFcn, defaultDriftFcn

from fpex0.InitialDistribution import InitialDistribution


def exampleSetup():
    """
    Generates example setup for FPEX0.

    ## Takes
    No input.

    ## Returns
    **setup** 
    <br> An FPEX0 setup configuration (`fpex0.setup.Setup`) example.
    <br>
    
    ### Comments
    Diffusion and drift (their value) is connected to the grid step size! <br>
    NOTE: So simulations and fits should be done with identical gridT.
    
    Initial and final diffusion are described by parameters!  
    --> Linear diffusion plus condition of a non-negative diffusion at end time
    
    Formula: 
    >       Diffusion = pD1 + t*(pD2-pD1)/betamax   
     
    where
    >       betamax --> maximum heat rate (maximum time in FP) ,
    >       pD1     --> initial diffusion at heat rate (time) 0
    >       pD2     --> final diffusion at heat rate (time) betamax

    NOTE:  Parameter bounds
    Parameter bounds should never be active in the solution. 

    Their only role is to limit search space and to ensure that no "invalid" valued are attained.
    """
    # Generate parameters object
    p0_FPdrift       = [  0.1  ,   0.1  ]    # linear drift
    p_lb_FPdrift     = [ -5.1  ,  -5.1  ]    # lower bounds
    p_ub_FPdrift     = [ 10.1  ,  10.1  ]    # upper bounds
    
    p0_FPdiffusion   = [  0.2  ,   0.1  ]    # linear diffusion
    p_lb_FPdiffusion = [  0.00 ,   0.00 ]    # lower bounds
    p_ub_FPdiffusion = [ 10.00 ,  10.00 ]    # upper bounds
    
    p0_iniDist       = [  2  ,  40  , 135 , 15.0 , 0.1000 ]  # Fraser-Sizuki parameters
    p_lb_iniDist     = [  1  ,  15  , 110 ,  0.1 , 1e-5   ]  # lower bounds
    p_ub_iniDist     = [  3  , 150  , 150 , 50.0 , 1.00   ]  # upper bounds
    
    parametersObj  = Parameters(    p0_FPdrift,   p0_FPdiffusion,   p0_iniDist,
                                    p_lb_FPdrift, p_lb_FPdiffusion, p_lb_iniDist,
                                    p_ub_FPdrift, p_ub_FPdiffusion, p_ub_iniDist  )

    # Generate grid 
    N        = 1001  # space resolution

    betamax  = 20   # maximum heat rate
    gridT    = np.linspace( 60,      160,  N )   # x-grid = temperatures
    gridTdot = np.linspace(  0,  betamax, 20 )   # t-grid = heating rates
    gridObj  = Grid(gridT, gridTdot)

    # set the FP functions and the initial distribution
    FPdriftFcn     = lambda t,p: FokkerPlanck.defaultDriftFcn(t,p)
    FPdiffusionFcn = lambda t,p: FokkerPlanck.defaultDiffusionFcn(t,p,betamax)

    # Fraser-Suzuki -> uses functions to describe the support
    FraserSuzuki_name   = "Fraser-Suzuki"
    FraserSuzuki_params = sympy.symbols("r,h,z,wr,sr")
    FraserSuzuki_expr   = sympy.sympify("h * exp( -log(r)/(log(sr)^2) * log(((x-z)*(sr^2-1))/(wr*sr) + 1)^2)")
    FraserSuzuki_zeropos = sympy.sympify("z - (wr*sr)/(sr^2-1)")
    FraserSuzuki_support = (-np.inf, FraserSuzuki_zeropos)
    
    FraserSuzuki = InitialDistribution(FraserSuzuki_expr, FraserSuzuki_support, FraserSuzuki_params, name=FraserSuzuki_name)
    IniDistFcn     = lambda x,p: FraserSuzuki.f(x,p)

    # Setup integrator (using defaults)
    integrationObj = Integration()

    # generate the setup object
    FPEX0setup = Setup(gridObj, parametersObj, integrationObj, FPdriftFcn, FPdiffusionFcn, IniDistFcn)

    return FPEX0setup


def importExampleMeasurements(FPEX0setup, gridskip: int):
    """ 
    Imports example measurement into a given setup FPEX0setup.
    <br>

    ## Takes
    **FPEX0setup**
    <br> An FPEX0 setup configuration (`fpex0.setup.Setup`).

    **gridskip**
    <br> Use grid with given stepsize (gridskip). 
    <br> For gridskip=2 for example, only every second grid point is used,
    whereas gridskip=1 uses whole grid.

    ## Returns
    **FPEX0setup**
    <br> The FPEX0setup passed in, now having imported example measurements in parameter FPEX0setup.Measurements.
    
    """
    for rate in ["0.60", "1.25", "2.50", "5.00", "10.00", "20.00"]:
        
        with resources.path('fpex0.example', f'ID407-rate_{rate}.json') as path:
            # using 'with' as resources.path() returns a context manager;
            # that context manager provides an os-path
            # --> path is now an os-path
            print(f"    {path}")
            with open(path) as file:
                # use again a context manager to let open() care about closing the file
                data = json.load(file)
                FPEX0setup.importMeasurements(data["T"], data["latentdata"], data["rate"], data["ID"], gridskip)

    return FPEX0setup