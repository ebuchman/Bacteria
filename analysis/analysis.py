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
        d = trajs[i][t+tau] - trajs[i][t]
        
        this_msd[tau] += d[0]**2 + d[1]**2
        
        this_msd_count[tau] += 1
        
    this_msd_count += this_msd_count==0
    this_msd = this_msd/this_msd_count
    msd += this_msd

  msd = msd/N
  
  return msd

def get_velocities(trajs, dt, skip, w):
  vs = np.zeros((trajs.shape[0], trajs.shape[1]-1, 2))

  for n in xrange(vs.shape[0]):
    for t in xrange(vs.shape[1]):
      vs[n][t][0] = min_sep(trajs[n][t+1][0] , trajs[n][t][0], w)
      vs[n][t][1] = min_sep(trajs[n][t+1][1] , trajs[n][t][1], w)


  return vs/(dt*skip)


def velocity_correlation(vs):
  TIME_POINTS = vs.shape[1]
  N = vs.shape[0]
  vc = np.zeros((TIME_POINTS))
  
  v2 = np.zeros((vs.shape[0], vs.shape[1]))

  for i in xrange(N):
    
    vc_count = np.zeros((TIME_POINTS))
    agent_vc = np.zeros((TIME_POINTS))
    
    for t in xrange(TIME_POINTS):
      # for each initial time, make a correlation curve.  run through all initial times, counting how many entries for each tau
      
      this_vc = np.zeros((TIME_POINTS))
      #v2 = np.mean(vs[:, t, 0]**2 + vs[:, t, 1]**2) # to normalize

      for tau in xrange(TIME_POINTS-t):
        vdot = np.dot(vs[i][t+tau], vs[i][t])
        this_vc[tau] += vdot
        vc_count[tau]+=1
      #this_vc = this_vc/v2
      agent_vc += this_vc
    agent_vc = agent_vc/vc_count
    vc += agent_vc

  vc = vc/N

  vc_norm = vc/np.mean(vs[:, :, 0]**2 + vs[:, :, 1]**2) # to normalize
  print np.mean(vs[:, :, 0]**2 + vs[:, :, 1]**2)
      
  return vc, vc_norm