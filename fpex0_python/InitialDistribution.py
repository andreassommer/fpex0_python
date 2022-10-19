import numpy as np
import sympy 
import numbers


class InitialDistribution:
   """Class holding initial distribution.

   ## Parameters

   **functionExpressionList** : list
   <br> list of piecewise function definitions

   **supportList** : list of (tuples or functions)
   <br> List of (open) support intervals, use numpy.inf and -numpy.inf for positive/negative infinity
   if a function is given, it will be evaluated with current parameters and must return a tuple for the support.

   **psymbols** : tuple of sympy symbols
   <br> Tuple of parameter set (used in all support intervals).

   **xsymbol** : sympy symbol
   <br> Symbol used as x variable (default: x).
         
   **name** : str
   <br> A handy name for the distribution.
   """

   def __init__(self, funcExprList, supportList, psymbols, xsymbol = sympy.symbols("x"), name = "unnamed"):
 
      # ensure the list variables are indeed lists (wrap singletons)
      if not(isinstance(funcExprList, list)):  funcExprList = [funcExprList]
      if not(isinstance( supportList, list)):   supportList = [supportList]
      
      # ensure we have as many functions as supports
      assert( len(funcExprList) == len(supportList) )

      # dimensions
      supportCount  = len(supportList)
      paramCount    = len(psymbols)

      # generate functions and derivatives w.r.t. p (i.e. dfdp)
      dfdp_sym = np.empty( (supportCount,paramCount) , dtype=object )
      dfdp_fun = np.empty( (supportCount,paramCount) , dtype=object )
      f_sym    = np.empty( supportCount              , dtype=object )
      f_fun    = np.empty( supportCount              , dtype=object )
      for (j , fExpr) in enumerate(funcExprList):
         f_sym[j] = funcExprList[j]                                 # store symbolic function  
         f_fun[j] = self.sympy2numpyXP(f_sym[j], xsymbol, psymbols) # transform to function
         for (i, p) in enumerate(psymbols):
            symfun = sympy.diff(fExpr, p)                                 # derivative w.r.t. p_i
            dfdp_sym[j,i] = symfun                                        # store symbolic function
            dfdp_fun[j,i] = self.sympy2numpyXP(symfun, xsymbol, psymbols) # transform to function

      # transform the support into support functions
      for (j, support) in enumerate(supportList):
         (l,r) = support               # extract left and right
         lfun = self.functionalizeSupportDescriptor(l, psymbols)
         rfun = self.functionalizeSupportDescriptor(r, psymbols)
         supportList[j] = (lfun,rfun)  # store support functions


      # store quantities (mostly for debugging)
      self.name     = name
      self.funcExprList  = funcExprList
      self.supportList   = supportList
      self.paramList     = psymbols
      self.xsymbol       = xsymbol
      self.supportCount  = supportCount
      self.paramCount    = paramCount
      self.f_fun    = f_fun
      self.f_sym    = f_sym
      self.dfdp_fun = dfdp_fun
      self.dfdp_sym = dfdp_sym

      # generate evaluator functions
      return None


   # some helpers 

   @staticmethod
   def getValuesWithinInterval(x, interval):
      """Returns values and indices of numpy array x that are inside given interval."""
      # NOTE: Boolean vs. integer indexing depends on size and count
      (l,r) = interval                # extract interval boundaries
      idx = np.where((x>l) & (x<r))   # returns a 1-element tuple 
      values = x[idx]                 # extract values
      return values, idx


   @staticmethod
   def evalSupport(support, p):
      """Evaluates support functions at given parameter vector."""
      (l,r) = support
      return (l(p), r(p))

   @staticmethod
   def sympy2numpyXP(fun, xsymbol, psymbols):
      """Lambdify specified function as a function of x and p."""
      xpsymbols = (xsymbol,) + tuple(psymbols)                # collect all symbols
      multiargfun = sympy.lambdify(xpsymbols, fun, 'numpy')   # function of x and individual parameters
      xpfun = lambda x,p : multiargfun(x, *p)                 # function that decollates vector p in multiple args
      return xpfun


   @staticmethod
   def sympy2numpyP(fun, symbols):
      """Lambdify specified function."""
      multiargfun = sympy.lambdify(tuple(symbols), fun, 'numpy') # function of individual parameters
      pfun = lambda p : multiargfun(*p)                          # function that decollated vector p in multiple args
      return pfun


   @staticmethod
   def functionalizeSupportDescriptor(thing, symbols):
      """Transforms a support descriptor into a callable function."""
      fun = None
      if isinstance(thing, sympy.Expr):      fun = InitialDistribution.sympy2numpyP(thing, symbols)
      if isinstance(thing, numbers.Number):  fun = lambda *ignore: thing
      if (fun == None): raise Exception("Invalid support descriptor", thing)
      return fun

   # @staticmethod
   # def compileFunction(fun):
   #    debugdir = "./codegen"
   #    cfun = autowrap(fun, tempdir=debugdir)
   #    return cfun


   def f(self, x, p):
      """Evaluates nominal function at x with parameter vector p.
      ## Takes
      **x**: numpy vector <br>
      **p**: numpy vector
      
      
      ## Returns:
      **f**: numpy vector 
      """
      fvals = np.zeros_like(x)
      for (j,support) in enumerate( self.supportList ):          # walk through the individual functions
         support = self.evalSupport(support, p)                  #   evaluate support functions
         (xx, idx) = self.getValuesWithinInterval( x , support ) #   get x-values in current support
         fvals[idx] = self.f_fun[j]( xx , p )                    #   evaluate function
      return fvals


   def dfdp(self, x, p):
      """Evaluates jacobian at x with parameter vector p.
      ## Takes
      **x**: numpy vector <br> 
      **p**: numpy vector


      ## Returns
      **J**: Numpy matrix dfdp, J[i][...] contains sensitivity w.r.t. p[i].
      """
      dfdpvals = np.zeros_like( x , shape = (len(p),len(x)) )        # preallocate
      for (j,support) in enumerate( self.supportList ):              # walk through supports
         support = self.evalSupport(support, p)                      #   determine support
         (x, idx) =  self.getValuesWithinInterval( x , support )     #   extract x values within support
         for i in range(0, self.paramCount):                         #   walk through derivative funcions
            dfdpvals[i][idx] = self.dfdp_fun[j,i](x,p)               #     evaluate derivative functions
      return dfdpvals

