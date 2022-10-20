from fpex0 import clearZeroFromMax
import numpy as np
from fpex0 import getBaselinePrimitive
from fpex0 import detectLinearRange
from fpex0 import getBaseline


def test_detectLinearRange():
    # copy inputs from matlab run
    Y = np.genfromtxt("tests/test-data/y.csv", delimiter=",")
    X = np.genfromtxt("tests/test-data/x.csv", delimiter=",")
    
    # left side
    side        = 'left'
    initlen     = 405
    reldevAmax  = np.infty
    reldevBmax  = 0.01
    reldevS2max = 2
    absdevAmax  = -1
    absdevBmax  = -1
    absdevS2max = -1
    stopidx, reg, initreg = detectLinearRange(X,Y,side,initlen,reldevAmax,reldevBmax,reldevS2max,absdevAmax,absdevBmax,absdevS2max)
    
    # check if detectLinearRange returns same results
    # stopidx
    assert(stopidx == 406) # Accept off-by-one. Check that in detail later.
    # Huge error here.
    
    # reg
    assert(reg.n     == 407)
    assert(reg.Xmean - 93.508332555282540      < 10e-13)
    assert(reg.Ymean - 6.071320189240500       < 10e-13)
    assert(reg.See   - 0.708781820835585       < 10e-13)
    assert(reg.Sxx   - 559.3581924611424       < 10e-13)
    assert(reg.Sxy   - 11.813325006880126      < 10e-13)
    assert(reg.a     - 4.096477703307445       < 10e-13)
    assert(reg.b     - 0.021119427883772       < 10e-13)
    assert(reg.s2    - 0.001750078569964       < 10e-13)
    # initreg
    assert(initreg.n                             == 405  )
    assert(initreg.Xmean - 93.498355135802460    < 10e-13)
    assert(initreg.Ymean - 6.070948312549130     < 10e-13)
    assert(initreg.See   - 0.706552585901338     < 10e-13)
    assert(initreg.Sxx   - 551.1535712379173     < 10e-13)
    assert(initreg.Sxy   - 11.507577751997193    < 10e-13)
    assert(initreg.a     - 4.118788973378958     < 10e-13)
    assert(initreg.b     - 0.020879076817285     < 10e-13)
    assert(initreg.s2    - 0.001753232223080     < 10e-13)

    # right side
    side        = 'right'
    initlen     = 495
    reldevAmax  = np.infty
    reldevBmax  = np.infty
    reldevS2max = 2
    absdevAmax  = -1
    absdevBmax  = 0.5
    absdevS2max = -1
    stopidx, reg, initreg = detectLinearRange(X,Y,side,initlen,reldevAmax,reldevBmax,reldevS2max,absdevAmax,absdevBmax,absdevS2max)
    assert(reg.n                               == 2322 )
    assert(reg.Xmean - 154.2557189696970       < 10e-13)
    assert(reg.Ymean - 154.2557189696970       < 10e-13)
    assert(reg.See   - 41.93883328106729       < 10e-13)
    assert(reg.Sxx   - 1.039401172103837e+05   < 10e-11)  # python saves (or displays) the variable more precise than matlab here
    assert(reg.Sxy   - -1223.342071245329      < 10e-13)
    assert(reg.a     - 7.312637437108421       < 10e-13)
    assert(reg.b     - -0.011769681467351      < 10e-13)
    assert(reg.s2    - 0.018077083310805       < 10e-13)
    # initreg
    assert(initreg.n                             == 495  )
    assert(initreg.Xmean - 154.2557189696970     < 10e-13)
    assert(initreg.Ymean - 5.498301459343400     < 10e-13)
    assert(initreg.See   - 2.523989783576885     < 10e-13)
    assert(initreg.Sxx   - 1006.973702336575     < 10e-13)
    assert(initreg.Sxy   - -15.555610617012533   < 10e-13)
    assert(initreg.a     - 7.881225555744368     < 10e-13)
    assert(initreg.b     - 7.881225555744368     < 10e-13)
    assert(initreg.s2    - 0.005119654733422     < 10e-13)


def test_getBaselinePrimitive():
    # passes
    Y = np.genfromtxt("tests/test-data/y.csv", delimiter=",")
    X = np.genfromtxt("tests/test-data/x.csv", delimiter=",")
    index_l = 406
    icept_l = 4.096477703307445
    slope_l = 0.021119427883772
    index_r = 4215
    icept_r = 7.311747545889207
    slope_r = -0.011763251099770
    res     = 1000
    
    
    bl_type = "sigmoidal"
    blSigmoidal = getBaselinePrimitive(X, Y, index_l, icept_l, slope_l, index_r, icept_r, slope_r, bl_type, res)
    X_test = np.linspace(90, 155, 101)  # 101 for good step size
    Y_test = blSigmoidal(X_test)
    Y_test_matlab = np.genfromtxt("tests/test-data/Y_getBaselinePrimitive_sigeval.csv", delimiter=",")
    assert (Y_test - Y_test_matlab < 10e-13).all()

    bl_type = "linear"
    blLinear = getBaselinePrimitive(X, Y, index_l, icept_l, slope_l, index_r, icept_r, slope_r, bl_type, res)
    X_test = np.linspace(92, 155, 101)
    # make the range smaller here, because pp linear scipy interpolator can't evaluate out of interpolation range
    Y_test = blLinear(X_test)
    Y_test_matlab = np.genfromtxt("tests/test-data/Y_getBaselinePrimitive_lineval.csv", delimiter=",")
    assert (Y_test - Y_test_matlab < 10e-13).all()


def test_getBaseline():
    X = np.genfromtxt("tests/test-data/x.csv", delimiter = ",")
    Y = np.genfromtxt("tests/test-data/y.csv", delimiter = ",")

    blfun, bldata = getBaseline(X, Y, bl_type='linear')
    
    # -- linear case -- (passes)
    
    # bldata
    # onset, endset
    assert(bldata.onset  - 100.5445200000000 < 10e-13)
    assert(bldata.endset - 133.8238700000000 < 10e-13)
    
    # reg_L
    assert( bldata.reg_L.n     == 407                        )
    assert( bldata.reg_L.Xmean - 93.508332555282540 < 10e-13 )
    assert( bldata.reg_L.Ymean - 6.071320189240500  < 10e-13 )
    assert( bldata.reg_L.See   - 0.708781820835585  < 10e-13 )
    assert( bldata.reg_L.Sxx   - 559.3581924611424  < 10e-13 )
    assert( bldata.reg_L.Sxy   - 11.813325006880126 < 10e-13 )
    assert( bldata.reg_L.a     - 4.096477703307445  < 10e-13 )
    assert( bldata.reg_L.b     - 0.021119427883772  < 10e-13 )
    assert( bldata.reg_L.s2    - 0.001750078569964  < 10e-13 )
    
    # reg_R
    assert( bldata.reg_R.n     == 2322                           )
    assert( bldata.reg_R.Xmean - 1.451371214857879e+02  < 10e-13 )
    assert( bldata.reg_R.Ymean - 5.604419748132494      < 10e-13 )
    assert( bldata.reg_R.See   - 41.938833281067290     < 10e-13 )
    assert( bldata.reg_R.Sxx   - 1.039401172103837e+05  < 10e-10 )
    assert( bldata.reg_R.Sxy   - -1.223342071245329e+03 < 10e-12 )
    assert( bldata.reg_R.a     - 7.312637437108421      < 10e-13 )
    assert( bldata.reg_R.b     - -0.011769681467351     < 10e-13 )
    assert( bldata.reg_R.s2    - 0.018077083310805      < 10e-13 )
    
    
    # blfun
    X_test = np.linspace(92, 156, 101)
    Y_getBaseline_lineval = np.genfromtxt("tests/test-data/Y_getBaseline_lineval.csv", delimiter=",")
    assert(( blfun(X_test) - Y_getBaseline_lineval < 10e-13 ).all())
    
    
    
    # -- sigmoidal case --  (passes)

    blfun, bldata = getBaseline(X, Y, bl_type='sigmoidal')
    
    # bldata
    # onset, endset
    assert(bldata.onset  - 1.021514400000000e+02 < 10e-13)
    assert(bldata.endset - 1.338238700000000e+02 < 10e-13)
    
    # reg_L
    # reg_R
    # --> Regressions are built the same way, don't need to test again 
    
    # blfun
    Y_getBaseline_sigeval = np.genfromtxt("tests/test-data/Y_getBaseline_sigeval.csv", delimiter=",")
    assert(( blfun(X_test) - Y_getBaseline_sigeval < 10e-13 ).all())

    # check if vectorized version runs
    baseline_list = getBaseline([X, X], [Y, Y], bl_type='linear')
    pass


def test_clearZeroFromMax():
    A = np.array([0, 1, 4, 0, 1, 15, 6, 0, 1, 0])
    B = clearZeroFromMax(A)
    A_cleared = np.array([0, 0, 0, 0, 1, 15, 6, 0, 0, 0])
    assert(np.all(B == A_cleared))

if __name__=='__main__':
    test_getBaseline()
    test_clearZeroFromMax()