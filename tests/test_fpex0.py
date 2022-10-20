import numpy as np

from fpex0 import simulate
from fpex0.example.exampleSetup import exampleSetup

def test_simulate():
    FPEX0setup = exampleSetup()

    # fraser suzuki parametrization
    r  =     2
    h  =    40
    z  = 135.0
    wr =    15
    sr =   0.1
    p_IC = [r, h, z, wr, sr]

    # FP drift and diffusion parametrization
    p_FPdrift     = [0.00000, 0.00000]
    p_FPdiffusion = [0.00000, 0.00000]

    # assemble p
    pvec = np.concatenate((p_FPdrift,  p_FPdiffusion,  p_IC))

    # NOTE: we could also use the initial parameter estimates from FPEX0setup
    # p = FPEX0setup.Parameters

    # simulate
    sol = simulate(FPEX0setup, pvec)

    # display solution
    # FPEX0_visualize(FPEX0setup, sol)

    # Plot.cla()
    # Plot.plot(sol.t, sol.y)
    # Plot.title(f"Solution")
    # Plot.show(block=False)

    return sol


if __name__=='__main__':
    sol = test_simulate()