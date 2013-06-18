import numpy as np
from openData import *


# traj shape is (num bacteria, time, xy)


def mean_squared_displacement(trajs):
  TIME_POINTS = trajs.shape[1]
  N = trajs.shape[0]
  msd = np.zeros((TIME_POINTS))
  
  for i in xrange(N):
    this_msd = np.zeros((TIME_POINTS))
    this_msd_count = np.zeros((TIME_POINTS)) #for averages over time
    
    for t in xrange(TIME_POINTS):
      for tau in xrange(TIME_POINTS-t):
        dx = trajs[i][t+tau][0] - trajs[i][t][0]
        dy = trajs[i][t+tau][1] - trajs[i][t][1]
        
        this_msd[tau] += dx**2 + dy**2
        
        this_msd_count[tau] += 1
    this_msd_count += this_msd_count==0
    this_msd = this_msd/this_msd_count
    msd += this_msd
  msd = msd/N
  
  return msd


def velocity_correlation(trajs):
  TIME_POINTS = trajs.shape[1]
  N = trajs.shape[0]
  vc = np.zeros((TIME_POINTS))

  
  
  
  return vc