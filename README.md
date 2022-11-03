# FPEX0 Python

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

The software comes with an example implemented in fpex0.example.exampleFit().
Running it will import example measurements, build an example setup and execute the algorithm.  
After about 20 steps it should give a solution near by (rounded):
>       p = [-0.9555,  0.03284, 0.2862, 3.4171, 2.5246, 43.0456, 131.8116, 3.5925, 0.1893]


## Extrapolating your own data

The heart of the package are the function fpex0.fit() and the class fpex0.Setup.
Setup holds your measurements and all problem-related configurations, e.g. your diffusion function and its
inital parameters. Then fit() uses your setup to fit the Fokker-Planck evolution to your measurements as
an optimization problem.
The example fit and used example functions should give a good understanding how the software is used.
Please also read about sympy symbolic functions if not familiar.


## Data processing

The functions described above assume **baseline corrected** data, so raw measurements must be processed.
The modules CP, baseline can do that for you.  
The processing consists of two parts: (1) calculating heat capacities, (2) detecting a baseline and
subtracting it.  
Both of it is done by addCP() plus some previous data preparation. As there is no code example, we
explain its usage:

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

Now you can modify the setup and extrapolate your data.

If you want to skip part (1) or (2), check for CP_DIN11357(), getBaseline(), subtractBaseline() or the
source code of addCP().


## About the implementation

This is a Python version of Andreas Sommer's matlab implementation, which can be found at
https://github.com/andreassommer/fpex0.

The Fokker-Planck equation is solved as an ODE via method of lines, using scipy
solve_ivp with BDF method as a default. That is basically a python version of matlab ode15s.  
The initial distribution, drift and diffusion are then fitted to the measurement data via an optimizer,
by default scipy least_squares (which is also currently the only option).  
Other optimizers and integrators can be implemented by the user, if compatible to the interplay of
fpex0.fit(), residual() and simulate(). However, the software is designed around the method of lines,
so using another method to solve Fokker-Planck will require significant adjustments.
