from fpex0 import fpex0
from fpex0.example.exampleSetup import exampleSetup, importExampleMeasurements

import time


def exampleFit():
    """
    Runs an example fit with:
    - An example fpex0 setup configuration, retrieved by `fpex0.example.exampleSetup.exampleSetup()` (check for details)
    - Example data, contained in the fpex0/example/ID407-rate_*.json files
    - A Fraser-Suzuki initial distribution
    - The scipy solve_ivp BDF integrator to solve the Fokker Planck equation
    - The scipy least_squares optimizer to fit the example measurements
    
    <br>
    No in- or output. Prints information and a solution parameter vector. <br>
    """
    # make example setup
    print("Creating an example setup.")
    FPEX0setup = exampleSetup()

    # import the example data
    print("Importing example data.")
    FPEX0setup = importExampleMeasurements(FPEX0setup, gridskip=2)

    # modify some configuration (as example)
    print("Setting integration options.")
    FPEX0setup.Integration.options["rtol"] = 1e-8
    FPEX0setup.Integration.options["atol"] = 1e-14

    # fitting
    print("\nSolving the fitting problem:")
    refTime = time.time()
    fit = fpex0.fit(FPEX0setup)
    duration = time.time() - refTime
    print("Fit complete.")
    print(f"Fit time: {duration}s.")
    print(f"Solution found: p = {fit.x}")

if __name__=='__main__':
    exampleFit()