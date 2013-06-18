import numpy as np
from matplotlib import pyplot as plt

from analysis import *

PATH = '/Users/BatBuddha/Programming/Research/Bacteria/code'


trajs, details = retrieve_data(PATH)
msd = mean_squared_displacement(trajs)
plt.plot(msd)
plt.show()