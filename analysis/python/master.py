import numpy as np
from matplotlib import pyplot as plt
import time

from analysis import *

PATH = '/Users/BatBuddha/Programming/Research/Bacteria/code/data/N1_T10000_G1_A0_dt100/'

trajs, details = retrieve_data(PATH)
vs = get_velocities(trajs, details['dt'], details['skip'], details['frame_lim'])

msd = msd_fast(trajs)**0.5
avgV = ensemble_avg_velocity(vs)
vc = velocity_correlation_np(trajs)



plt.subplot(311)
plt.plot(msd)
plt.xscale('log')
plt.yscale('log')
plt.subplot(312)
plt.plot(avgV)
plt.subplot(313)
plt.plot(vc)

plt.figure()
plt.hist(avgV)
plt.show()
quit()


