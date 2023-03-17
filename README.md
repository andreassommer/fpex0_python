# FPEX0 Python

Authors: Michael Strik and Andreas Sommer at the Interdisciplinary Center for Scientific Computing (IWR), Heidelberg University.

Documentation: https://andreassommer.github.io/fpex0_python/

This package gives a Python implementation of the FPEX0 method  
for data-driven de-smearing of DSC signals presented in the paper

Sommer, Andreas; Hohenauer, Wolfgang; Barz, Tilman:  
Data-driven de-smearing of DSC signals.  
J Therm Anal Calorim (2022).  
https://doi.org/10.1007/s10973-022-11258-y


A Matlab version of FPEX0 is available at
https://github.com/andreassommer/fpex0



## Installing the package

The fpex0 package can be installed via pip:
```
pip install fpex0
```



## Running an example

The software comes with an example implemented in `fpex0.example.exampleFit()`, that can be run by:
```python
import fpex0
fpex0.example.exampleFit()
```
It will import measurements, build a setup and execute the algorithm.  
After about 20 steps it should give a solution near by:
>       p = [-0.9555,  0.03284, 0.2862, 3.4171, 2.5246, 43.0456, 131.8116, 3.5925, 0.1893]


## Extrapolating your own data

The heart of the package are the function `fpex0.fit()` and the class `fpex0.Setup`.
`Setup` holds your measurements and all problem-related configurations, e.g. your initial distribution and 
its inital parameters. Then `fit()` uses your setup to fit the Fokker-Planck evolution to the measurements
as an optimization problem.  
Reading the source code of `exampleFit()` and `exampleSetup()` should give a good understanding how the 
software can be used. We also recommend reading about sympy symbolic functions if not familiar.


## Data processing

The functions described above assume **baseline corrected** data, so raw measurements must be processed.
The submodules `CP`, `baseline` can do that for you.  
Processing consists of two parts: 
* calculating heat capacities,
* detecting a baseline and subtracting it.

Both of it is done by `addCP()`, it will also do some previous data preparation.  
As there is no code example, we will explain its usage:

1. Create a DSC_Data object and load measurements
```python
dsc_data = DSC_Data()
dsc_data.T = T
dsc_data.dsc = dsc
dsc_data.rate = rate
```

2. Process
```python
dsc_data = addCP(dsc_data)
```

3. Create fpex0 setup and import your data
```python
FPEX0setup = Setup(gridObj, parametersObj, integrationObj, FPdriftFcn, FPdiffusionFcn, IniDistFcn)
FPEX0setup.importDSCobj(dsc_data)
```

Now you can modify the setup and finally extrapolate your data via
```python
fit = fpex0.fit(FPEX0setup)
```

Note that both `addCP()` and `importDSCobj()` can also take lists of `DSC_Data` objects:
```python
dsc_data = addCP([dsc_data1, dsc_data2, dsc_data3])
# now dsc_data holds a list
FPEX0setup.importDSCobj(dsc_data)
```

If you want to skip one of the steps, check for
* `CP_DIN11357()` and `CP_sapphire_DIN11357()` (heat capacities),
* `subtractBaseline()` and `getBaseline()` (baseline).

These sets of functions can execute the raw single steps.


## About the implementation

This is a Python version of Andreas Sommer's matlab implementation, which can be found at
https://github.com/andreassommer/fpex0.

The Fokker-Planck equation is solved as an ODE via method of lines, using scipy's
solve_ivp with BDF method as a default. This is basically a python version of matlab ode15s.  
The initial distribution, drift and diffusion are then fitted to the measurement data via an optimizer,
by default scipy's least_squares (which is also currently the only option).  
Other optimizers and integrators can be implemented by the user, if compatible to the interplay of
`fpex0.fit()`, `residual()` and `simulate()`. However, the software is designed around the method of lines,
so using another method to solve Fokker-Planck will require significant adjustments.
