# this is where I can mess with sparse matrices to see how they are manipulated

import numpy as np
from scipy import sparse
import time
import sys
from scipy.sparse import random

data1 = np.array([[2,0,0,0,5],
                  [0,0,0,3,0],
                  [0,0,7,0,0],
                  [4,5,6,7,0],
                  [0,8,0,2,0]])
data2 = np.array([[2, 0, 0, 0, 5],
                  [0, 0, 0, 3, 0],
                  [0, 0, 7, 0, 0],
                  [4, 5, 6, 7, 0],
                  [0, 8, 0, 2, 0]])
data3 = np.array([[np.inf, 0, 0, 0, 5],
                  [0, 0, 0, 3, 0],
                  [0, 0, 7, 0, 0],
                  [4, 5, 6, 7, 0],
                  [0, 8, 0, np.inf, 0]])











