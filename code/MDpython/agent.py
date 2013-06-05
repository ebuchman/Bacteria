import numpy as np
import matplotlib.pyplot as plt
from pylab import rand

BALL_R = 0.2
NUM_BALLS = 3
FRAME_LIM = 10

gamma = 1.

class Agent():
    def __init__(self, cm, v0, ball_r, theta = None, num_balls = NUM_BALLS, frame_lim = FRAME_LIM, dt= 0.01):
        
        self.r = ball_r
        self.N = num_balls
        self.cm = cm
        self.v0 = v0
        self.F_self = 1.
        self.dt = dt

        self.has_neighbours = False

        self.frame_lim = frame_lim
        
        if theta == None:
            self.th = np.random.rand()*2*np.pi                
        else:
            self.th = theta
            
        self.omega = 0
        self.v = 0 #self.v0
        
        
        self.t = 0  #time since last cell division
        
        self.balls = self.compute_rod() #initialize hard rod



        
        self.last_F = 0
        self.last_tau = 0
        
    
    def compute_rod(self):
        balls = []
        for i in xrange(self.N):
            x = (self.cm[0] - (self.N - abs(self.N%2 - 1) - 2*i)*self.r*np.cos(self.th))%self.frame_lim
            y = (self.cm[1] - (self.N - abs(self.N%2 - 1) - 2*i)*self.r*np.sin(self.th))%self.frame_lim
            
            balls.append([x, y])
        return balls
            
    
    
    def step(self, F, tau):
        #velocity verlet
    
        
        F += self.F_self - gamma*self.v
        tau += (np.random.rand()-0.5)*2*2*np.pi - gamma*self.omega
        
        #update v and omega
        dv = 0.5*(F + self.last_F)*self.dt
        domega = 0.5*(tau + self.last_tau)*self.dt

        self.v += dv
        self.omega += domega
                
        #update cm and th
        dx = self.v*np.cos(self.th)*self.dt + 0.5*F*np.cos(self.th)*self.dt**2
        dy = self.v*np.sin(self.th)*self.dt + 0.5*F*np.sin(self.th)*self.dt**2
        
        dth = self.omega*self.dt + 0.5*tau*self.dt**2
        
        self.cm[0] = (self.cm[0] + dx)%self.frame_lim
        self.cm[1] = (self.cm[1] + dy)%self.frame_lim
        self.th = (self.th + dth)%(2*np.pi)
        
        self.balls = self.compute_rod()
        
        
        self.last_F = F
        self.last_tau = tau
        
        self.t += 1
        
        
        

