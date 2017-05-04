import sys
from SALib.sample import saltelli
from scipy import stats
import numpy as np

def test():
    m = np.mean(np.arange(10))
    sys.stdout.write('testing SAVIO!! %0.2f' % (m))
    
if __name__ == '__main__':
    test()