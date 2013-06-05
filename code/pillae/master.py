from agent import Agent
from flock import Flock

culture = Flock(num_agents = 5, run_time = 100000, show_every = 10)
culture.run()
