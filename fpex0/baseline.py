import numpy as np
from scipy import interpolate
from scipy import integrate
from copy import copy
from fpex0.linreg import LinearRegression


class BaselineData:
    """
    Stores information about the baselevel.
    
    ## Parameters
    **reg_L**
    <br> Left linear part.
    
    **reg_R**
    <br> Right linear part.
    
    **onset** 
    <br> Onset of phase transition.
    
    **endset**
    <br> Endset of the phase transition.
    """
    def __init__(self, reg_L=None, reg_R=None, onset=None, endset=None):
        self.reg_L  = reg_L
        self.reg_R  = reg_R
        self.endset = endset
        self.onset  = onset

class BaselineDetectionSettings:
    """
    Stores (default) baseline detection settings to be passed to `fpex0.baseline.detectLinearRange()`.
    
    ## Parameters
    Naming: 
    - R,L = right or left
    - abs, rel = absolute or relative
    - dev = deviation
    - initfraction = initial fraction of data to interpolate.

    **LabsdevA**, **LabsdevB**, **LabsdevS2**
    
    **LreldevA**, **LreldevB**, **LreldevS2**
    
    **Linitfraction**

    **RabsdevA**, **RabsdevB**, **RabsdevS2**
    
    **RreldevA**, **RreldevB**, **RreldevS2**
    
    **Rinitfraction**
    """
    def __init__(self):
        self.LabsdevA      = -1
        self.LabsdevB      = -1
        self.LabsdevS2     = -1
        self.LreldevA      = np.infty
        self.LreldevB      = 0.01
        self.LreldevS2     = 2.00
        self.Linitfraction = 0.1
        
        self.RabsdevA      = -1
        self.RabsdevB      = 0.05
        self.RabsdevS2     = -1
        self.RreldevA      = np.infty
        self.RreldevB      = np.infty
        self.RreldevS2     = 2.00
        self.Rinitfraction = 0.2


def detectLinearRange(X, Y, side, initlen, reldevAmax, reldevBmax, reldevS2max, absdevAmax, absdevBmax, absdevS2max):
    """
    Detect the stop position of the linear range.

    ## Takes
    
    **X**: x-values
        
    **Y**: y-values
           
    **side**: 'L' or 'R' for "from left" or "from right"
        
    **initlen**
    <br> Length (in samples) to calculate the initial standard deviation.
         
    **reldevA** 
    <br> Acceptable relative deviation of initial y-axis intercept.
     
    **reldevB**
    <br> Acceptable relative deviation of initial slope.
     
    **reldevS2**
    <br> Acceptable relative deviation of initial squared standard deviation.
     
    **absdevA**
    <br> Acceptable absolute deviation of initial y-axis intercept.
     
    **absdevB**
    <br> Acceptable absolute deviation of initial.
     
    **absdevS2**
    <br> Acceptable absolute deviation of initial squared standard deviation.


    ## Returns
    **stopidx** 
    <br> Position in X, where the linear range is estimated to be left.
     
    **reg**
    <br> Final regression structure as returned by fpex0.linreg.LinearRegression.linReg.
     
    **initreg**
    <br> Initial regression structure as returned by fpex0.linreg.LinearRegression.linReg.
 
    ### Comments
     * Set the `reldev` to inf to disable
     
     * Set the `absdev` to  -1 to disable
     
     * The `reldev` is only be checked if value exceeds `absdev`
    
    """
    
    
    # choose loop-range
    if side.lower() == "left":
        start = initlen
        stop  = len(X)
        step  = 1
        # initial regression
        initreg = LinearRegression()
        initreg.linReg(X[:start], Y[:start])
    elif side.lower() == "right":
        start = len(X)-1 - initlen  # start at position left from initial range
        stop  = -1                        # exclusive index 
        step  = -1
        # initial regression
        initreg = LinearRegression()
        initreg.linReg(X[start+1:], Y[start+1:])


    reg = copy(initreg)
    index = None

    for index in range(start, stop, step):
        # update regression
        reg.linReg(X[index], Y[index])  

        absdevA  = abs(reg.a - initreg.a    )   
        reldevA  = abs(absdevA  / initreg.a )
        absdevB  = abs(reg.b - initreg.b    )   
        reldevB  = abs(absdevB  / initreg.b )
        absdevS2 = abs(reg.s2 - initreg.s2  )  
        reldevS2 = abs(absdevS2 / initreg.s2)

        # check deviation
        if (absdevA  > absdevAmax ) and (reldevA  > reldevAmax ):
            break
        if (absdevB  > absdevBmax ) and (reldevB  > reldevBmax ):
            break
        if (absdevS2 > absdevS2max) and (reldevS2 > reldevS2max):
            break
    stopidx = index

    return stopidx, reg, initreg


def getBaselinePrimitive(X, Y, index_l, icept_l, slope_l, index_r, icept_r, slope_r, bl_type, res=1000):
    """
    Retrieves the baselevel function for specified data.
    
    ## Takes
      
    **X**: x-values (temperatures, time, ...)
     
    **Y**: y-values (voltage, heat flux, ....)
     
    **index_l**
    <br> Index where left linear part is left ("onset")
     
    **icept_l**
    <br> y-intercept of left linear part.
     
    **slope_l**
    <br> Slope of left linear part.
     
    **index_r**
    <br> Index where right linear part is left ("onset").
     
    **icept_r**
    <br> y-intercept of right linear part.
     
    **slope_r** 
    <br> Slope of right linear part.
     
    **bl_type** 
    <br> What type should the baselevel function have? "linear" or "sigmoidal".
     
    **res**
    <br> Number of support points for sigmoidal (default: 100).
    

    ## Returns
      
    **blfun**
    <br> Function handle of baselevel function
     
    """

    # calculate values at "onset" and "offset"
    yval_l = icept_l + X[index_l] * slope_l
    yval_r = icept_r + X[index_r] * slope_r
    yval_delta = yval_r - yval_l
    
    # calculate linear base level function
    lin_bl_xvals = [        X[0]        , X[index_l], X[index_r],          X[-1]         ]
    lin_bl_yvals = [icept_l+X[0]*slope_l,   yval_l  ,   yval_r  ,  icept_r+X[-1]*slope_r ]
    lin_bl_fun  = interpolate.interp1d(lin_bl_xvals, lin_bl_yvals, 'linear')  # piecewise linear interpolation
    
    
    if bl_type.lower() == 'linear':
        return lin_bl_fun
    
    elif bl_type.lower() in ['sigmoid', 'sigmoidal']:
        
        # sigmoidal interpolation as in Fig. 3 of DIN 51007

        # substract baseline and integrate peak part "in between"
        Xmid = X[index_l:index_r+1]
        Ymid = Y[index_l:index_r+1] - lin_bl_fun(Xmid)
        Ymid = np.array([y if y >= 0 else 0 for y in Ymid])  # set negatives to zero
        cumarea = integrate.cumtrapz(Ymid, Xmid)  # cumulative integral (area)
        sigmoidal = np.concatenate(([0], cumarea)) / max(cumarea) * yval_delta + yval_l  # add a zero for consistency with matlab code

        # interpolate integral at support points (res = #intervals)
        sig_nodes = np.linspace(Xmid[0], Xmid[-1], res+1) # Here should be an error in the corresponding matlab code. It produces 1001 nodes for res=1000
        sig_interp = interpolate.interp1d(Xmid, sigmoidal, 'linear')  # this is now a function object
        sig_nodevals = sig_interp(sig_nodes)

        # generate baseline function (piecewise cubic hermite interpolant)
        sig_x = np.concatenate( ( [X[0]], [X[index_l-1]], sig_nodes, [X[index_r+1]],  [X[-1]]) )
        sig_y = np.concatenate( ( [lin_bl_fun( X[0] )], [lin_bl_fun( X[index_l-1] )], sig_nodevals, [lin_bl_fun( X[index_r+1] )],  [lin_bl_fun( X[-1] )] ) )
        blfun = interpolate.PchipInterpolator(sig_x, sig_y)  # this is now a function object
        
        return blfun



def getBaseline(X, Y, bl_type='linear', blds=BaselineDetectionSettings()):
    """
    Retrieves the baselevel function for specified DSC data.
    
    ## Takes
    **X**: 
    <br> Numpy vector of x-values (e.g. temperatures) or list (!) of vectors.
     
    **Y**: 
    <br> Numpy vector of y-values (e.g. cp-values) or list (!) of vectors.
     
    **bl_type**: 'linear' or 'sigmoidal'
    (default: 'linear')
     
    **blds** 
    <br> BaselineDetectionSettings-object containing the BaseLine Detection Setting
    (default: Default BLDS object).
    

    ## Returns
    **blfun**
    <br> Function handle of baselevel function.
     
    **bldata**
    <br> Structure containing information about the baseline.

    OR if list have been passed:

    **baseline_list**:
    <br> List of tuples (blfun, bldata).
    """
    
    
    if type(X) is list:
        baseline_list = [getBaseline(X[i], Y[i], bl_type, blds) for i in range(len(X))]
        return baseline_list
    
    print("Detecting linear ranges")
    peakPos = np.argmax(Y)
    initlen_L = int( np.floor(blds.Linitfraction * peakPos) )
    initlen_R = int( np.floor(blds.Rinitfraction * (len(X)-peakPos)) )
    idx_L, reg_L, _ = detectLinearRange(X,Y,'left' ,initlen_L,blds.LreldevA, blds.LreldevB, blds.LreldevS2, blds.LabsdevA, blds.LabsdevB, blds.LabsdevS2)
    idx_R, reg_R, _ = detectLinearRange(X,Y,'right',initlen_R,blds.RreldevA, blds.RreldevB, blds.RreldevS2, blds.RabsdevA, blds.RabsdevB, blds.RabsdevS2)
    print('Done.')

    # get baselevel function
    blfun = getBaselinePrimitive(X, Y, idx_L, reg_L.a, reg_L.b, idx_R, reg_R.a, reg_R.b, bl_type)
    
    # small plausibility check:  slope of linear parts should be small;
    maxslope = 0.1
    if (abs(reg_L.b) > maxslope) or (abs(reg_R.b) > maxslope):
       print("\n")
       print(f"Slope of linear part is large: left: {reg_L.b}, right: {reg_R.b}. There's probably something wrong. Using proposed baseline instead!")
       if abs(reg_L.b) > maxslope:
           reg_L.b = 0
           reg_L.a = Y[0]
           idx_L = 1
       if (abs(reg_R.b) > maxslope): 
           reg_R.b = 0
           reg_R.a = Y[-1]
           idx_R = len[Y]-2
       blfun = getBaselinePrimitive(X, Y, idx_L, reg_L.a, reg_L.b, idx_R, reg_R.a, reg_R.b, bl_type)

    
    # build data
    
    # classic onset/endset estimation: where the linear parts are left (poor estimation)
    # bldata.onset  = X[idx_L]
    # bldata.endset = X[idx_R]
    
    # better onset/endset estimation: point X where the data (X,Y) first falls below baseline (seen from peak maximum)
    # with fallback using the ind_L or idx_R
    bloffset = 0.02
    temp = Y[0:peakPos+1] - blfun(X[0:peakPos+1])
    counter = None
    # get index of last value below bloffset
    for k in range(len(temp)):
        if temp[k] < bloffset:
            counter = k
    idxOnset = counter
    
    temp = Y[peakPos:] - blfun(X[peakPos:])
    counter = None
    # get index of first value below bloffset
    for k in range(len(temp)):
        if temp[k] < bloffset:
            counter = k
            break
    idxEndset = counter + peakPos - 1
    
    if idxOnset is None:
        idxOnset  = idx_L
    if idxEndset is None:
        idxEndset = idx_R
    
    # save data
    bldata = BaselineData()  # save all together at the end?
    bldata.reg_L  = reg_L  # regression information left
    bldata.reg_R  = reg_R  # regression information right    
    bldata.onset  = X[idxOnset]
    bldata.endset = X[idxEndset]

    return blfun, bldata


def subtractBaseline(X, Yin, blfun, onset=None, endset=None, clearzero=True, nonnegative=True):
    """
    Subtracts the baseline from given points.
    
    ## Takes   
    **X**: x values (e.g. vector temperatures)
    
    **Y**: y values or function (e.g. vector of cp values, or function cp(T))
    
    **blfun**
    <br> Function handle to baseline function.
    
    **clearzero**
    <br> Flag indicating to clear zeros (see DSC204_clearZeroFromMax) (default: True).
    
    **nonnegative**
    <br> Flag indicating to ensure nonnegativity                      (default: True).
    
    **onset**
    <br> Onset value (zero values are put below/left of this x value) [optional].
    
    **endset**
    <br> Endset value (zero values are put above/right of this x value) [optional].
    

    ## Returns   
    **Yvals** 
    <br> Processed y values.

    **Yfun**
    <br> Interpolator of processed values.
    """

    if onset is None:
        onset = min(X)
    if endset is None:
        endset = max(X)
    
    assert( len(X) == len(Yin) )
    Yvals = Yin

    # subtract baseline from Y data
    Yvals = Yvals - blfun(X)

    # make zeros outside the interval [onset, endset]
    idx = [i for i in range(len(X)) if (X[i] < onset or X[i] > endset)]
    Yvals[idx] = 0

    # ensure nonnegativity
    if nonnegative:
        idx = [i for i in range(len(Yvals)) if Yvals[i] < 0]
        Yvals[idx] = 0

    # clear zero
    if clearzero:
       Yvals = clearZeroFromMax(Yvals)

    
    # for the interpolation function, add some zero space left and right
    addlen = 5
    XX = np.concatenate( ( X[0] + np.arange(-addlen,0), X       , X[-1] + np.arange(1,addlen+1) ) )
    YY = np.concatenate( ( np.zeros(addlen)             ,  Yvals  , np.zeros(addlen)            ) )

    # build function from values
    Yfun = interpolate.PchipInterpolator(XX,YY)

    return Yvals, Yfun


def clearZeroFromMax(data):
    """ 
    Locates the maximum value in datavec, searches the first occurances
    of a zero value to the left and to the right from the peak position
    and deletes the noise before and after these positions.
    
    ## Takes 
    **data**: Data array to be cleared
    
    ## Returns
    **data**: Cleared data array
    """
    # Find index/position of maximal value
    maxidx = np.argmax(data)
    
    # Find the positions of the first zeros, seen from the peak value
    clearToIdx   = max( np.where(data[:maxidx+1]==0)[0] )
    clearFromIdx = min( np.where(data[maxidx:]  ==0)[0] ) + maxidx
    
    # delete noise
    if not (clearToIdx is None):
        data[:clearToIdx] = 0
    if not (clearFromIdx is None):
        data[clearFromIdx:] = 0

    return data
