import numpy as np
import copy

from fpex0 import DSC_Data
from fpex0 import addCP

def test_addCP():

    sample = DSC_Data()
    sample.T = np.genfromtxt("tests/test-data/T_ID16-407.csv")
    sample.dsc = np.genfromtxt("tests/test-data/dsc_ID16-407.csv")
    sample.mass = 10.47
    sample.rate = 0.6

    reference = DSC_Data()
    reference.T = np.genfromtxt("tests/test-data/T_REF-ID16-398s.csv")
    reference.dsc = np.genfromtxt("tests/test-data/dsc_REF-ID16-398s.csv")
    reference.mass = 84.73
    reference.rate = 10
    sample_copy = copy.deepcopy(sample)
    sample = addCP(sample, reference)
    cp = sample.cp
    

    # import and check matching values
    cp_values           = np.genfromtxt("tests/test-data/cp_values_ID16-407.csv")
    cp_T                = np.genfromtxt("tests/test-data/cp_T_ID16-407.csv")
    cp_latentdata       = np.genfromtxt("tests/test-data/cp_latentdata_ID16-407.csv")
    cp_latentfun_eval   = np.genfromtxt("tests/test-data/cp_latentfun_eval100-120_ID16-407.csv", delimiter=",")
    assert( min ( cp.values - cp_values < 10e-12 ) )
    assert( min ( cp.T - cp_T < 10e-12 ) )
    assert( min ( cp.latentdata - cp_latentdata < 10e-12 ) )
    assert( min ( cp.latentfun(np.arange(100,121)) - cp_latentfun_eval < 10e-12 ) )

    # check if vectorized version runs
    sample_copy2 = copy.deepcopy(sample_copy)
    samples = addCP([sample_copy, sample_copy2], reference)
    
if __name__=='__main__':
    test_addCP()
