from copy import copy

import numpy as np
from scipy.integrate import solve_ivp
from scipy import interpolate

from .fokker_planck import FokkerPlanck


class DSC_Data:
    """
    Class for dsc data. Stores unprocessed and processed data.
    
    ## Parameters
    **mass**
    <br> Sample mass.

    **rate**
    <br> Heating rate.
    
    **T**
    <br> Temperatures.
    
    **dsc**
    <br> Corresponding raw DSC measurements.
    
    **bldata** 
    <br> Information about the baselevel.
    
    **blfun** 
    <br> Baselevel function as retrieved by `fpex0.baseline.getBaseline`.
    
    **cp** 
    <br> Calculated heat capacities as HeatCapacitiy member.

    **ID**
    <br> Name/identifier of experiment.
    """
    def __init__(self, T=None, dsc=None, mass=None, rate=None, bldata=None, blfun=None, cp=None, ID=None):
        self.mass      = mass     
        self.rate      = rate     
        self.T         = T        
        self.dsc       = dsc      
        self.bldata    = bldata   
        self.blfun     = blfun
        self.cp        = cp
        self.ID        = ID



class Setup:
    """
    FPEX0 configuration.

    ## Parameters : `fpex0.setup.Grid`
    **Grid**
    <br> A discretization grid for Fokker-Planck.
    
    **Parameters**
    <br> A collection of parameters and boudns, 
    instance of `fpex0.setup.Parameters`.
    
    **Integration** : `fpex0.setup.Integration`
    <br> Represents an integrator to solve our discretization of Fokker-Planck.

    **FPdriftFcn** 
    <br> Function object of Fokker-Planck drift.
    
    **FPdiffusionFcn**
    <br> Function object of Fokker-Planck diffusion.
    
    **IniDistFcn**
    <br> Function object of inital distribution.
    """
    def __init__(self, Grid, Parameters : dict, Integration, FPdriftFcn, FPdiffusionFcn, IniDistFcn):
        self.Grid           = Grid
        self.Parameters     = Parameters
        self.Integration    = Integration
        self.FPdriftFcn     = FPdriftFcn
        self.FPdiffusionFcn = FPdiffusionFcn
        self.IniDistFcn     = IniDistFcn
        
        self.Measurements   = Measurements()


    def make_rhsFcn(self, p_FPdrift, p_FPdiffusion):
        """
        Generator for the ODEs right-hand-side.
        """
        h = self.Grid.h
        FPdriftFcn = self.FPdriftFcn
        FPdiffusionFcn = self.FPdiffusionFcn

        rhsFcn = lambda t,u: FokkerPlanck.FokkerPlanckODE(t, u, self.Grid.h, self.FPdriftFcn, p_FPdrift, self.FPdiffusionFcn, p_FPdiffusion)
        return rhsFcn


    def make_jacFcn(self, p_FPdrift, p_FPdiffusion):
        """
        Generator for the ODEs jacobian.
        """
        fokker_planck = FokkerPlanck()
        jacFcn = lambda t,u: fokker_planck.FokkerPlanckODE_dfdu(t, u, self.Grid.h, self.FPdriftFcn, p_FPdrift, self.FPdiffusionFcn, p_FPdiffusion)
        return jacFcn


    def importMeasurements(self, T, values, heatrate, ID=None, gridskip=1):
        """
        Imports measurements into FPEX0setup (`fpex0.setup.Measurements`).
        An interpolation to the integration grid is performed and the interpolated values are stored.
        
        ## Takes 
        **FPEX0setup**
        <br> FPEX0_class_setup object.
        
        **ID** 
        <br> Name/identifier of experiment.
        
        **heatrate**
        <br> Heating rate.
        
        **T**
        <br> Temperature grid / measurement positions.
        
        **values**
        <br> Measurement values (e.g. cp-values); 
        where values[k] denotes the measurement value at T[k].
        
        **gridskip**
        <br> Use grid with given stepsize (gridskip). 
        <br> For gridskip=2 for example, only every second grid point is used,
        whereas gridskip=1 uses the whole grid.
        

        ## Returns
        No return value.
        """
        
        # ensure that the rate (=time) lies inside the gridTdot
        if heatrate < self.Grid.gridTdot[0] or heatrate > self.Grid.gridTdot[-1]:
           print(f'importMeasurements: Specified rate {heatrate} outside grid bounds [{self.Grid.gridTdot[0]}, {self.Grid.gridTdot[-1]}] !')
        
        # ensure we have a single column
        rawT    = np.reshape(T     , -1)  # passing -1 creates a 1d array
        values  = np.reshape(values, -1)

        # interpolate values to gridT respecting gridskip
        rawValues = values
        temperatures = self.Grid.gridT[::gridskip]
        interp = interpolate.PchipInterpolator(rawT, rawValues, extrapolate=False)
        values = interp(temperatures)
        values = np.where(np.isnan(values), 0, values)   # replace NaN by 0
        
        # add measurements
        self.Measurements.appendMeasurement(rawT, temperatures, rawValues, values, heatrate, ID)


    def importDSCobject(self, DSC_Data, gridskip=1):
        # T, values, heatrate, ID=None, gridskip=1

        if type(DSC_Data) is list:
            for i in range(len(DSC_Data)):
                self.importDSCobject(DSC_Data[i], gridskip)
            return

        T = DSC_Data.T
        # use cp data if available, raw dsc otherwise
        if DSC_Data.cp is not None:
            values = DSC_Data.cp
        else: 
            values = DSC_Data.dsc
            print("WARNING: Using unprocessed dsc data. Make sure it is baseline-corrected!")
        heatrate = DSC_Data.rate
        ID = DSC_Data.ID

        self.importMeasurements(T, values, heatrate, ID, gridskip)


class Parameters:
    """
    ## Parameters

    **p0**
    <br> Full initial parameter vector.

    **p_lb**
    <br> Full parameter lower bound vector.
    
    **p_ub**
    <br> Full parameter upper bound vector.

    **p0_FPdrift**
    <br> Initial drift parameter vector.

    **p0_FPdiffusion**
    <br> Initial diffusion parameter vector.
    
    **p0_iniDist**
    <br> Initial paramater vector for initial distribution.
    
    **p_lb_FPdrift**
    <br> Lower bound vector for drift parameters.

    **p_lb_FPdiffusion**
    <br> ...
    
    **p_lb_iniDist**
    <br> ...
    
    **p_ub_FPdrift**
    <br> Upper bound vector for drift parameters.
    
    **p_ub_FPdiffusion**
    <br> ...
    
    **p_ub_iniDist**
    <br> ...

    **idxFPdrift**
    <br> Indices of drift parameters in p0 (resp. p_lb, p_ub).

    **idxFPdiffusion**
    <br> Same thing, diffusion parameters. 
    
    **idxFPall**
    <br> Same thing, all Fokker-Planck parameters.
    
    **idxIniDist**
    <br> Same thing, initial distribution parameters.
    """
    def __init__(self,  p0_FPdrift,        p0_FPdiffusion,        p0_iniDist,
                        p_lb_FPdrift=None, p_lb_FPdiffusion=None, p_lb_iniDist=None,
                        p_ub_FPdrift=None, p_ub_FPdiffusion=None, p_ub_iniDist=None ):

        # ensure matching lengths
        assert(p_lb_FPdrift is None or len(p0_FPdrift) == len(p_lb_FPdrift) 
           and p_ub_FPdrift is None or len(p0_FPdrift) == len(p_ub_FPdrift))

        assert(p_lb_FPdiffusion is None or len(p0_FPdiffusion) == len(p_lb_FPdiffusion) 
           and p_ub_FPdiffusion is None or len(p0_FPdiffusion) == len(p_ub_FPdiffusion))

        assert(p_lb_iniDist is None or len(p0_iniDist) == len(p_lb_iniDist) 
           and p_ub_iniDist is None or len(p0_iniDist) == len(p_ub_iniDist))

        # initial parameters
        self.p0_FPdrift         = p0_FPdrift
        self.p0_FPdiffusion     = p0_FPdiffusion
        self.p0_iniDist         = p0_iniDist

        # lower parameter bounds
        # set to vectors of -np.inf if not given
        if p_lb_FPdrift is None:
            self.p_lb_FPdrift = np.full(len(p0_FPdrift), -np.inf)
        else:
            self.p_lb_FPdrift = p_lb_FPdrift

        if p_lb_FPdiffusion is None:
            self.p_lb_FPdiffusion = np.full(len(p0_FPdiffusion), -np.inf)
        else:
            self.p_lb_FPdiffusion = p_lb_FPdiffusion

        if p_lb_iniDist is None:
            self.p_lb_iniDist = np.full(len(p0_iniDist), -np.inf)
        else:
            self.p_lb_iniDist = p_lb_iniDist

        # upper parameter bounds
        # set to vectors of np.inf if not given
        if p_ub_FPdrift is None:
            self.p_ub_FPdrift = np.full(len(p0_FPdrift), np.inf)
        else:
            self.p_ub_FPdrift = p_ub_FPdrift

        if p_ub_FPdiffusion is None:
            self.p_ub_FPdiffusion = np.full(len(p0_FPdiffusion), np.inf)
        else:
            self.p_ub_FPdiffusion   = p_ub_FPdiffusion

        if p_ub_iniDist is None:
            self.p_ub_iniDist = np.full(len(p0_iniDist), np.inf)
        else:
            self.p_ub_iniDist       = p_ub_iniDist

        # put initials, lower and upper bounds in one vector each
        self.p0 = np.concatenate((p0_FPdrift, p0_FPdiffusion, p0_iniDist))
        self.p_lb = np.concatenate((p_lb_FPdrift, p_lb_FPdiffusion, p_lb_iniDist))
        self.p_ub = np.concatenate((p_ub_FPdrift, p_ub_FPdiffusion, p_ub_iniDist))
        
        # indices
        self.idxFPdrift      = np.arange(0, len(p0_FPdrift) )
        self.idxFPdiffusion  = np.arange(0, len(p0_FPdiffusion) ) + len(p0_FPdrift)
        self.idxFPall        = np.arange(0, len(p0_FPdrift)+len(p0_FPdiffusion) )
        self.idxIniDist      = np.arange(0, len(p0_iniDist) ) + ( len(p0_FPdrift) + len(p0_FPdiffusion) )


    # extractor functions for external parameter vector (alike self.p0)
    def extract_p_all(self, pvec):
        return pvec

    def extract_p_FPdrift(self, pvec):
        return pvec[self.idxFPdrift]

    def extract_p_FPdiffusion(self, pvec):
        return pvec[self.idxFPdiffusion]

    def extract_p_FPall(self, pvec):
        return pvec[self.idxFPall]

    def extract_p_iniDist(self, pvec):
        return pvec[self.idxIniDist]


    # setter functions (TODO: Set parameters in the vectors.)
    def set_p0(self, p0):
        self.p0 = p0

    def set_p_lb(self, p_lb):
        assert(len(p_lb) == len(self.p0))
        self.p_lb = p_lb
    
    def set_p_ub(self, p_ub):
        assert(len(p_ub) == len(self.p0))
        self.p_ub = p_ub



class Grid:
    """
    FPEX0 discretization(?) grid.

    ## Parameters
    **gridT**
    <br> x-grid = temperatures
    <br> This is our method of lines discretization grid.
    
    **gridTdot**
    <br> t-grid = heating rates
    
    **N** 
    <br> Number of grid points.
    
    **h** 
    <br> Grid size.
    """
    def __init__(self, gridT, gridTdot):
        self.gridT = gridT
        self.gridTdot = gridTdot
        self.N = len(gridT)
        self.h = ( gridT[-1] - gridT[0] ) / self.N



class Measurements:
    """
    Stores measurement data in list containers.
    
    ## Parameters
     **rawT** 
     <br> Measurement temperatures.
     
     **temperatures** 
     <br> Selected temperatures from temperature grid.
     
     **rawValues** 
     <br> Raw masurement values.
     
     **values** 
     <br> rawValues interpolated to grid.
     
     **rates** 
     <br> Heat rates.
     
     **ID** 
     <br> Experiment name/identifier.
    """
    def __init__(self, rawT=[], temperatures=[], rawValues=[], values=[], rates=[], ID=[]):
        self.rawT           = rawT
        self.temperatures   = temperatures
        self.rawValues      = rawValues
        self.rates          = rates
        self.ID             = ID
        self.values         = values

    def appendMeasurement(self, new_T, new_temperatures, new_rawValues, new_values, new_rate, new_ID=None):
        """
        Appends given measurements to the measurement parameters. 
        """
        self.rawT.append(new_T)
        self.temperatures.append(new_temperatures)
        self.rawValues.append(new_rawValues)
        self.values.append(new_values)
        self.rates.append(new_rate)
        self.ID.append(new_ID)
        


class Integration:
    """
    Represents FPEX0 integrator.

    ## Parameters
    **method**
    <br> Method to use by scipy solve_ivp. Default is BDF.
    
    **options**
    <br> Integration options. Currently not used. 
    
    **integrator**
    <br> Integrator for ODEs. Should have the signature integrator(FPrhs, t0tf, u0, method).

    **NN, A, B**: Internal variables used by FokkerPlanckODE_dfdu.
    """
    # (TODO) Set default tolerances in Integration class (?)
    def __init__(self, integrator=solve_ivp, method="BDF", options={}):
        self.method     = method
        self.options    = options
        self.integrator = integrator
        self.NN         = None  # variables for FokkerPlanckODE_dfdu()
        self.A          = None
        self.B          = None

    def updateOptions(self, newOptions):
        """ Returns the integrators options. """
        options = copy(self.options)
        for key in newOptions:
            options[key] = newOptions[key]
        self.options = options
        return options
    

    def updateJacobian(self, jacobianFcn):
        """ Returns an option dictionary with updated jacobian. """
        options = copy(self.options)
        options["jac"] = jacobianFcn
        self.options = options
        return options
        

