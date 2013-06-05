import numpy as np
from pylab import rand

BALL_R = 0.2
NUM_BALLS = 3
FRAME_LIM = 10

NUM_PILLAE = 8
PILLI_SPAN = np.pi/3

k = 2.
Ff = 5.

gamma = 1.

class Pillus():
    def __init__(self):
        self.L0 = 0 #initial length
        self.L = 0 #retracted length
        self.F = 1.
        self.th = 0
        self.P = 10. #motor power


class Agent():
    def __init__(self, cm, ball_r, theta = None, num_balls = NUM_BALLS, num_pillae = NUM_PILLAE, pilli_span = PILLI_SPAN, frame_lim = FRAME_LIM, dt= 0.01):
        
        self.r = ball_r
        self.N = num_balls
        self.cm = cm
        self.F_self = 1.
        self.dt = dt

        self.num_pillae = num_pillae   
        self.pilli_span = pilli_span 
        self.pilli_length_mean = 2 # lengths are gaussian around this
        self.pilli_length_std = 0.1
         
        # list of pilli.  Each pillus has a length and force 
        self.pilli = []
        for i in range(num_pillae):
            self.pilli.append(Pillus())

        self.has_neighbours = False

        self.frame_lim = frame_lim
        
        if theta == None:
            self.th = np.random.rand()*2*np.pi                
        else:
            self.th = theta

        
        self.t = 0  #time since last cell division
        
        self.balls = self.compute_rod() #initialize hard rod



        self.v = 0
        self.omega = 0
        self.last_F = 0
        self.last_tau = 0
        
    
    def compute_rod(self):
        balls = []
        for i in xrange(self.N):
            x = (self.cm[0] - (self.N - abs(self.N%2 - 1) - 2*i)*self.r*np.cos(self.th))%self.frame_lim
            y = (self.cm[1] - (self.N - abs(self.N%2 - 1) - 2*i)*self.r*np.sin(self.th))%self.frame_lim
            
            balls.append([x, y])
        return balls
            
    def extend_pillus(self, pillus):
        pillus.L0 = np.random.normal(self.pilli_length_mean, self.pilli_length_std)
        pillus.L = pillus.L0 - 0.1
        pillus.th = np.random.uniform(-self.pilli_span/2., self.pilli_span/2.)
    
    def compute_pilli_forces(self):
        F_r, F_n, Tau  = 0, 0, 0
        
        for p in self.pilli:
            x = p.L0-p.L
            
            if x > 0:
                p.F = k*x
                
                if p.F < Ff:
                    #retract
                    if not p.F == 0: 
                        p.L -= (p.P/p.F)*self.dt
                else:
                    p.F = Ff
                    F_r += p.F*np.cos(p.th)
                    F_n += p.F*np.sin(p.th)
                    if not p.F == 0:
                        p.L -= (p.P/p.F)*self.dt

                Tau += p.L*np.sin(p.th)*p.F

        
        return F_r, F_n, Tau
                
    
    def step(self, F, tau):
        #velocity verlet
        for p in self.pilli:
            #print p.L
            if p.L <= 0:
                if np.random.rand() > 0.9:
                    #print 'extended'
                    self.extend_pillus(p)

                
        pilli_forces = self.compute_pilli_forces()
        #print pilli_forces
        F += pilli_forces[0] - gamma*self.v
        F_n = pilli_forces[1]
        tau += pilli_forces[2] - gamma*self.omega
        
        #update v and omega
        dv = 0.5*(F + self.last_F)*self.dt
        domega = 0.5*(tau + self.last_tau)*self.dt

        self.v += dv
        self.omega += domega
                
        #update cm and th
        dx = self.v*np.cos(self.th)*self.dt + 0.5*F*np.cos(self.th)*self.dt**2
        dy = self.v*np.sin(self.th)*self.dt + 0.5*F*np.sin(self.th)*self.dt**2
        
        dth = self.omega*self.dt + 0.5*tau*self.dt**2
        
        for p in self.pilli:
            p.th -= dth
        
        self.cm[0] = (self.cm[0] + dx)%self.frame_lim
        self.cm[1] = (self.cm[1] + dy)%self.frame_lim
        self.th = (self.th + dth)%(2*np.pi)
        
        self.balls = self.compute_rod()
        
        
        self.last_F = F
        self.last_tau = tau
        
        self.t += 1
        
        #print self.pilli[0].L
        

