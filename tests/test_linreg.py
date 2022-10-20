import numpy as np
from fpex0 import LinearRegression

# Test module: pytest

def test_LinReg():
    
    X = np.array([1, 2, 3, 4])
    Y = np.array([1, 3, 4, 8])
    reg = LinearRegression()
    
    # initial regression
    reg.linReg(X, Y)
    # values taken from a run of the matlab code
    assert(reg.n     - 4    < 10e-13)
    assert(reg.Xmean - 2.5  < 10e-13)
    assert(reg.Ymean - 4    < 10e-13)
    assert(reg.See   - 1.8  < 10e-13)
    assert(reg.Sxx   - 5    < 10e-13)
    assert(reg.Sxy   - 11   < 10e-13)
    assert(reg.a     - -1.5 < 10e-13)
    assert(reg.b     - 2.2  < 10e-13)
    assert(reg.s2    - 0.9  < 10e-13)
    
    # scalar update
    reg.linReg(5, 10)
    assert(reg.n     - 5    < 10e-13)
    assert(reg.Xmean - 3    < 10e-13)
    assert(reg.Ymean - 5.2  < 10e-13)
    assert(reg.See   - 1.9  < 10e-13)
    assert(reg.Sxx   - 10   < 10e-13)
    assert(reg.Sxy   - 23   < 10e-13)
    assert(reg.a     - -1.7 < 10e-13)
    assert(reg.b     - 2.3  < 10e-13)
    assert(reg.s2    - 0.63333333333333 < 10e-13)

    # vector update
    X = np.array([6, 7])
    Y = np.array([20, 23])
    reg.linReg(X, Y)
    
    assert(reg.n     - 7    < 10e-13)
    assert(reg.Xmean - 4    < 10e-13)
    assert(reg.Ymean - 9.857142857142858  < 10e-13)
    assert(reg.See   - 37.571428571428590  < 10e-13)
    assert(reg.Sxx   - 28   < 10e-13)
    assert(reg.Sxy   - 106   < 10e-13)
    assert(reg.a     - -5.285714285714285 < 10e-13)
    assert(reg.b     - 3.785714285714286  < 10e-13)
    assert(reg.s2    - 7.514285714285718 < 10e-13)
    
        
