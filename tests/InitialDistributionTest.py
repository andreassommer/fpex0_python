from fpex0 import InitialDistribution
import sympy
import numpy as np
import time
import matplotlib.pyplot as Plot


## Prepare some initial distributions

# Gaussian
Gaussian_name    = "Gaussian"
Gaussian_x       = sympy.symbols("x")
Gaussian_params  = sympy.symbols("mu,sigma")
Gaussian_expr    = sympy.sympify("1 / (sigma*sqrt(2*pi))*exp(-0.5*(x-mu)^2/(sigma^2))")
Gaussian_support = (-np.inf, np.inf)
Gaussian_pvals   = [0, 1]
Gaussian = InitialDistribution(Gaussian_expr, Gaussian_support, Gaussian_params, Gaussian_x, name=Gaussian_name)

# cutted Gaussian
cGaussian_name    = "cuttedGaussian"
cGaussian_support = (-1, 1)
cGaussian = InitialDistribution(Gaussian_expr, cGaussian_support, Gaussian_params, Gaussian_x, name=cGaussian_name)

# Fraser-Suzuki -> uses functions to describe the support
FraserSuzuki_name   = "Fraser-Suzuki"
FraserSuzuki_params = sympy.symbols("r,h,z,wr,sr")
FraserSuzuki_expr   = sympy.sympify("h * exp( -log(r)/(log(sr)^2) * log(((x-z)*(sr^2-1))/(wr*sr) + 1)^2)")
FraserSuzuki_zeropos = sympy.sympify("z - (wr*sr)/(sr^2-1)")
FraserSuzuki_support = (-np.inf, FraserSuzuki_zeropos)
FraserSuzuki_pvals   = [2, 40, 135, 15, 0.1]  # r h z wr sr
FraserSuzuki = InitialDistribution(FraserSuzuki_expr, FraserSuzuki_support, FraserSuzuki_params, name=FraserSuzuki_name)



def testDistribution(distribution, pvals, psyms, l, r, c):
   print("Processing ", distribution.name)
   xvals = np.linspace(l,r,c)
   pvals = np.array(pvals)
   # ===== Evaluate nominal f
   sTime = time.time()
   fvals = distribution.f(xvals, pvals)
   duration = time.time() - sTime
   print(np.c_[xvals,fvals])
   print("p", pvals)
   print(distribution.name, "- Nominal took", duration, "seconds.")
   Plot.cla()
   Plot.plot(xvals, fvals)
   Plot.title("%s - Nominal f - Took: %s" % (distribution.name, duration))
   Plot.show(block=False)
   input("Press Enter to continue...")
   # ===== Now the jacobian dfdp
   sTime = time.time()
   dfdp = distribution.dfdp(xvals, pvals)
   duration = time.time() - sTime
   print(distribution.name, "- Jacobian took", duration, "seconds.")
   for (i,psym) in enumerate(psyms):
      Plot.cla()
      Plot.plot(xvals, dfdp[i])
      Plot.title("%s - df/d%s = dfdp_%d - Took: %ss" % (distribution.name, psyms[i], i, duration))
      Plot.show(block=False)
      input("Press Enter to continue...")
   
   


testDistribution(Gaussian, Gaussian_pvals, Gaussian_params, l=-2, r=+2, c=141)
testDistribution(cGaussian, Gaussian_pvals, Gaussian_params, l=-2, r=+2, c=141)
testDistribution(FraserSuzuki, FraserSuzuki_pvals, FraserSuzuki_params, l=-10, r=150, c=2000)

print("Done.")