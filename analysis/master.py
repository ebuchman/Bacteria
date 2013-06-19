import numpy as np
from matplotlib import pyplot as plt

from analysis import *

PATH = '/Users/BatBuddha/Programming/Research/Bacteria/code'


trajs, details = retrieve_data(PATH)

vs = get_velocities(trajs, details['dt'], details['skip'], details['frame_lim'])

'''
msd = mean_squared_displacement(trajs)
plt.plot(msd)
plt.show()
quit()
'''

vc, vc_norm = velocity_correlation(vs)
plt.subplot(211)
plt.plot(vc)
plt.subplot(212)
plt.plot(vc_norm)
plt.show()