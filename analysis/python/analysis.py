import numpy as np
from openData import *

#########################################################################################
# Fast Versions of Analysis routines (using numpy correlation functions)
#########################################################################################

# traj shape is (num bacteria, time, xy)

#########################################################################################
# Mean Squared Displacement Function
# splitting up the squared difference gives three terms.  Two are sums over squares of positions.  The other is a correlation
def msd_fast(trajs):
  T = trajs.shape[1]
  N = trajs.shape[0]

  trajs2 = trajs**2
  
  msd = np.zeros((T))
  
  for n in xrange(N):
    r2 = np.zeros((T))
    rtau2 = np.zeros((T))
    
    # compute sums over squares of positions for r(t) and r(t+tau)
    for tau in xrange(T):
      r2[tau] = np.sum(trajs2[n][:(T-tau), 0] + trajs2[n][:(T -tau), 1])
      rtau2[tau] = np.sum(trajs2[n][tau:, 0] + trajs2[n][tau:, 1])

    # compute auto correlation
    corx = np.correlate(trajs[n][:, 0], trajs[n][:, 0], mode='full')[T-1:]
    cory = np.correlate(trajs[n][:, 1], trajs[n][:, 1], mode='full')[T-1:]
    cor = corx + cory
    
    msd += (rtau2 - 2*cor + r2)  / np.arange(T, 0, -1)

  msd = msd/N

  return msd

#########################################################################################
# Velocity Correlation Function
# normalized by product of v(t) and v(t+tau)
def velocity_correlation_np(vs):
  T = vs.shape[1]
  N = vs.shape[0]
  from matplotlib import pyplot as plt
  
  vs2 = vs**2
  #print vs2.shape
  magV = (vs2[:, :, 0] + vs2[:, :, 1])**0.5

  vs = vs/magV[:, :, None]

  vc = np.zeros((T))
  for i in xrange(N):

    this_vc_x = np.correlate(vs[i][:, 0], vs[i][:, 0], mode='full')[T-1:]
    this_vc_y = np.correlate(vs[i][:, 1], vs[i][:, 1], mode='full')[T-1:]
    
    this_vc = this_vc_x + this_vc_y
    vc += this_vc
  
  vc = vc / np.arange(T, 0, -1) # normalize for number of contributions per tau 
  vc = vc/N # normalize for number of bacteria

  return vc

#########################################################################################

#########################################################################################
# Velocity related functions
#########################################################################################

def get_velocities(trajs, dt, skip, w):
  vs = np.zeros((trajs.shape[0], trajs.shape[1]-1, 2))

  for n in xrange(vs.shape[0]):
    for t in xrange(vs.shape[1]):
      vs[n][t][0] = min_sep(trajs[n][t+1][0] , trajs[n][t][0], w)
      vs[n][t][1] = min_sep(trajs[n][t+1][1] , trajs[n][t][1], w)

  vs = vs/(dt*skip)

  return vs

def ensemble_avg_velocity(vs):
  vs = np.mean(vs, 0)
  vs = (vs[:, 0]**2 + vs[:, 1]**2)**0.5
 
  return vs
  
#########################################################################################
# Slow Versions of Analysis routines (using nested python loops)
# these can be sped up with numba, but the numpy routines above are exceedingly faster
#########################################################################################

def msd_0(trajs):
  T= trajs.shape[1]
  N = trajs.shape[0]

  msd = np.zeros((T))
  for n in xrange(N):
    for tau in xrange(T):
      d = trajs[n][tau] - trajs[n][0]
      d = d[0]**2 + d[1]**2
      msd[tau] += d
  msd = msd/N

  return msd


def msd_py(trajs):
  T = trajs.shape[1]
  N = trajs.shape[0]
  
  msd = np.zeros((T)) #velocity correlation as function of tau
  
  for n in xrange(N):    
    for tau in xrange(T):
      # for each tau: add contributions from all t0
      t_msd = 0
      
      for t in xrange(T-tau):
        
        d = trajs[n][t+tau] - trajs[n][t]
        d = d[0]**2 + d[1]**2
        t_msd += d

      t_msd = t_msd/(T-tau)
      msd[tau] += t_msd
  
  msd = msd/N
  
  return msd

def velocity_correlation_py(vs):
  TIME_POINTS = vs.shape[1]
  N = vs.shape[0]
  
  vc = np.zeros((TIME_POINTS)) #velocity correlation as function of tau

  for n in xrange(N):
    for tau in xrange(TIME_POINTS):
      # for each tau: add contributions from all t0
      t_vc = 0
      for t in xrange(TIME_POINTS-tau):
        cor = np.dot(vs[n][t+tau], vs[n][t])
        magVt = (vs[n][t+tau, 0]**2 + vs[n][t+tau, 1]**2)**0.5
        magV = (vs[n][t, 0]**2 + vs[n][t, 1]**2)**0.5
        t_vc += cor/magVt/magV
      t_vc = t_vc/(TIME_POINTS-tau)
      vc[tau] += t_vc

  vc = vc/N

  return vc
