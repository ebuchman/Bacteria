from agent import Agent
from flock import Flock
import cPickle as pickle

import time

t0 = time.clock()
culture = Flock(num_agents = 15, run_time = 1000, show_every = 10, animation = 'off')
trajectories, details = culture.run()
t1 = time.clock()

print 'python took ', t1-t0


f = open('trajectories.pkl', 'wb')
pickle.dump([trajectories, details], f)
f.close()
