import numpy as np
from matplotlib import pyplot as plt
import time

from analysis import *

PATH = '/Users/BatBuddha/Programming/Research/Bacteria/code/data/N1_T10000/'


trajs, details = retrieve_data(PATH)

vs = get_velocities(trajs, details['dt'], details['skip'], details['frame_lim'])
norm2 = np.mean(np.mean(vs[:, :, 0]**2 + vs[:, :, 1]**2, 1), 0)

msd = msd_0(trajs)
avgV = ensemble_avg_velocity(vs)
vc = velocity_correlation_np(vs)/norm2

plt.subplot(311)
plt.plot(msd)
plt.subplot(312)
plt.plot(avgV)
plt.subplot(313)
plt.plot(vc)
plt.show()
quit()


