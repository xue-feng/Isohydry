import sys
from SALib.sample import saltelli
from scipy import stats
import numpy as np
from params_soil import soil_dict
import matplotlib.pyplot as plt
import multiprocessing as mp
import cPickle as pickle

def test():
    m = np.mean(np.arange(10))
    soil_type = 'loamy_sand'
    Ps = soil_dict[soil_type]['Ps']
    sys.stdout.write('testing SAVIO!! The mean of array is %0.2f, Psat_soil is %.3e' % (m, Ps))
    with open('savio_testing.pickle', 'wb') as handle:
        pickle.dump(m, handle)
    pickle.dump()
    
if __name__ == '__main__':
    test()