from scipy import interpolate
from numpy.polynomial import Polynomial
import numpy as np

from fpex0.baseline import getBaseline
from fpex0.baseline import subtractBaseline



class HeatCapacity():
   """
   Stores cp, usually calculated by fpex0.CP.addCP().
   
   ## Parameters
   **T**: Temperatures
   
   **values**: Corresponding heat capacities
   
   **fun**: Interpolant (of T, values) (function handle).
   
   **latentdata**
   <br> Values with subtracted baseline.
   
   **latentfun** 
   <br> Interpolant of latentdata as (function handle).
   """
   def __init__(self, T=None, values=None, fun=None, latentdata=None, latentfun=None):
      self.T            = T
      self.values       = values
      self.fun          = fun
      self.latentdata   = latentdata
      self.latentfun    = latentfun

def CP_DIN11357(T, mS, dscS, cpR, mR, dscR, dsc0=0):
   """

   Calculates the (apparent) specific heat capacity.
   
   It applies the "heat flow calibration" method: A known reference cp (of sapphire) 
   is rescaled using the mass-ratio and signal ratio of reference and sample.
    
   >     cpS(T) = cpR(T) * mR/mS * (dscS(T)-dsc0(T)) / (dscR(T)-dsc0(T))

   ### Reference
   DIN EN ISO 11357-4:2014 <br>
   Plastics - Differential scanning calorimetry <br>
   Determination of specific heat capacity
    

   ## Takes

   **T**: Vector of temperatures to evaluate the cp-value
   
   **mS**: Mass of sample
   
   **dscS**: DSC signal of sample (microvolts)
   
   **cpR**: cp-values of reference (cf. fpex0.CP_sapphire_DIN11357())
   
   **mR**: Mass of reference

   **dscR**: DSC signal of reference (microvolts)
   
   **dsc0**: DSC signal with two empty crucibles
   

   ## Returns
   **cpS**: cp-values of sample at specified temperatures.
   
   ### Comments 
     * If vectors are given, they must coincide in size and meaning, <br>
              i.e.  dscS[k], dscR[k], dsc0[k], cpR[k] all correspond to temperature T[k]
     * If dscS and dscR are already baseline-corrected, set dsc0 to 0
   
   """
   # calculate cp
   cpS = cpR * mR/mS * (dscS-dsc0) / (dscR-dsc0)

   return cpS


def addCP(DSCsample, DSCreference, DSC0=0, Tmin=55, Tmax=160, bltype='linear'):
   """
   Adds cp values to DSC data structure (setup.DSC_Data).
   Calculation of cp is done with CP_DIN11357.
   
   ## Takes
   **DSCsample** 
   <br> fpex0.setup.DSC_Data object or list of objects, sample - (e.g. pcm)

   **DSCreference**
   <br> fpex0.setup.DSC_Data object, reference (e.g. saphire)
   
   **DSC0** 
   <br> DSC signal with two empty crucibles.                 [default: 0]
   <br> Only needed if signals are not zero-corrected yet.
   
   **Tmin**
   <br> Lower temperature bound.                             [default:  55 degC]
   
   **Tmax**
   <br> Upper temperature bound.                             [default: 160 degC]
   
   
   ## Returns
   **DSCsample**
   <br> fpex0.setup.DSC_Data object or list of objects with additional field cp.
   """
   
   # (TODO) Make these controllable
   scaleByMass = True  # multiply signaly by mass (signals are often normed to mass 1)
   scaleToRate = True  # scale signals to heatrate (higher rates induce higher signals, see below)
   TrangeWarn  = True  # warn if Tmin or Tmax exceed sample's or reference's temperature range
   
   # make function applicable for struct arrays
   if type(DSCsample) is list:
      DSCsample = [addCP(DSCsample[i], DSCreference, DSC0, Tmin, Tmax, bltype) for i in range(len(DSCsample))]
      return np.array(DSCsample)
   
   # retrieve the heating rates
   betaS = DSCsample.rate
   betaR = DSCreference.rate # possibly a vector
   
   # reference rates not unique?
   if type(betaR) is np.ndarray and len(betaR) != len(set(betaR)):
      print(f'DSC204_addCP: Reference heat rates not unique! rate vector = {betaR}\n')
      raise ValueError('Reference heat rates not unique!')

   print(f'DSC204_addCP: Processing ID={DSCsample.ID}')
   # Variable naming: AB
   #    A:    T -> Temperature   m -> mass   dsc -> dsc data
   #    B:    S -> Sample    R -> Reference  

   # quick accessors
   TR   = DSCreference.T
   TS   = DSCsample.T
   dscS = DSCsample.dsc
   dscR = DSCreference.dsc
   
   # calculate minima and maxima of sample and reference temperatures
   minTS = min(TS)
   maxTS = max(TS)
   minTR = min(TR)
   maxTR = max(TR)
   
   # issue a warning if Tmax or Tmin exceeds reference or sample temperature range
   if ( Tmin < max(minTS,minTR)  or  Tmax > min(maxTS,maxTR) ):
      if TrangeWarn:
         print('WARNING: Specified Tmin or Tmax exceeds sample or reference temperature range. Will be adjusted!')
    
    
   # determine Tmin and Tmax
   Tmin = max( [minTS, minTR, Tmin] )
   Tmax = min( [maxTS, maxTR, Tmax] )
    
   # restrict temperatures and align everything at the temperature information of the sample measurement
   # note: this is the temperature of the empty reference crucible, not of the sample itself
   # get the signals of the reference corresponding to the temperatures of the sample.
   # we use linear interpolation and extrapolation here
   idxR = [i for i in range(len(TR)) if Tmin <= TR[i] and TR[i] <= Tmax]
   idxS = [i for i in range(len(TS)) if Tmin <= TS[i] and TS[i] <= Tmax]
   dscR = dscR[idxR]
   dscS = dscS[idxS]
   TR   = TR[idxR]
   TS   = TS[idxS]
   dscR_interp = interpolate.interp1d(TR, dscR, kind='linear', fill_value='extrapolate')
   dscR = dscR_interp(TS)
                              
   # masses
   mS = DSCsample.mass     # mass of sample
   mR = DSCreference.mass  # mass of reference
    
   # measurements are normalized to uV/mg, so we recover the original signal by multiplying with mass
   if scaleByMass:
      dscS = mS * dscS
      dscR = mR * dscR
    
   # from carefully looking at the measurement data, we see that the voltage signal is proportional
   # to the heating rate, with proportionality constant approximately 1.
   # so we normalize both the sample and the reference signal to a heating rate of 1.0 K/min.
   # NOTE: this does also not interfere if betaR==betaS
   if scaleToRate:
      dscS = dscS / betaS
      dscR = dscR / betaR
    
   # now retrieve the reference cp values of saphire (unit: degC)
   cpR = CP_sapphire_DIN11357(TS, 'degC')
    
   # calculate cp of sample
   cpS = CP_DIN11357(TS, mS, dscS, cpR, mR, dscR, DSC0)
    
   # store cp values and associated temperatures
   cp = HeatCapacity()
   cp.values = cpS
   cp.T      = TS
   
   # store the piecewise polynomial with pchip interpolation; python checks bounds (Tmin, Tmax) for us
   cp.fun     = interpolate.PchipInterpolator(TS, cpS)
    
   # get the baseline
   blfun, bldata = getBaseline(TS, cpS, bltype)
    
   # build the latent cp function  subtractBaseline(X,  Yin, blfun, onset,        endset       , clearzero, nonnegative)
   latentCPvals, latentCPfun =     subtractBaseline(TS, cpS, blfun, bldata.onset, bldata.endset, False    , True       )
   cp.latentdata  = latentCPvals
   cp.latentfun   = latentCPfun
   
   # store it in DSC data
   DSCsample.cp = cp
   DSCsample.bldata = bldata
   DSCsample.blfun = blfun
   return DSCsample


def CP_sapphire_DIN11357(T, unit='degC'):
   """
   Delivers specific heat capacity of saphire for specified temperature according to DIN EN ISO 11357-4:2014-10.
   
   ## Takes
   **T**: Temperature <br>
   **unit**: one of "degC" or "K"    [default: degC]
   
   ## Returns
   **cp**
   <br> Literature value of sapphire as functor.
   

   [!] *Note*: This approximation is only valid in the interval 100K < T < 1200K.
   """

   # coefficients in J/(g*K)
   A = [
       1.12705,
       0.23260,
      -0.21704,
       0.26410,
      -0.23778,
      -0.10023,
       0.15393,
       0.54579,
      -0.47824,
      -0.37623,
       0.34407
   ]

   # linear transformations of temperature
   if unit.lower() in ['degc', 'c', 'celsius']:
      x = (T - 376.85) / 550
   elif unit.lower() in ['k', 'kelvin']:
      x = (T - 650) / 550
   else:
      raise ValueError(f"Unknown unit: {unit}")

   # build and evaluate polynomial
   cp_poly = Polynomial(A)
   cp = cp_poly(x)

   return cp
