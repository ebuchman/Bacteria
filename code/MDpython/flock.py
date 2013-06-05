import numpy as np
import matplotlib.pyplot as plt

from agent import Agent
from drawing import set_canvas, draw_

FRAME_LIM = 10
BALL_R = 0.2
NUM_BALLS = 3
DIVISION_LENGTH = NUM_BALLS*2
V0 = 1.

E = 1.
L = BALL_R*2

GROW_TIME = 20000


dt = 0.01

class Flock():
    def __init__(self, num_agents, run_time, show_every, animation):
        self.animation = animation
        
        self.N = num_agents
        self.run_time = run_time
        self.show_every = show_every
                
        self.frame_lim = FRAME_LIM
        
        self.time = 0
        
        self.dt = dt
        
        self.grow_time = GROW_TIME
        self.division_length = DIVISION_LENGTH
        
        self.agents = []
        for i in xrange(self.N):
            cm = np.random.rand(2)*FRAME_LIM
            new_agent = Agent(cm = cm, v0 = V0, ball_r = BALL_R, num_balls = NUM_BALLS, frame_lim = FRAME_LIM)
            self.agents.append(new_agent)
             
    def compute_forces(self, i, j, F_r, Tau):
        #loop over all ball-ball interactions
        for a in xrange(self.agents[i].N):
            for b in xrange(self.agents[j].N):
                dx = self.agents[i].balls[a][0] - self.agents[j].balls[b][0]
                dy = self.agents[i].balls[a][1] - self.agents[j].balls[b][1]
                r2 = dx**2 + dy**2
                
                if r2 < (2**(1./6)*L)**2:
                    F_piece = (48*E)*((L**12)*(r2**(-6)) - 0.5*(L**6)*(r2**(-3)))/r2
                    
                    F_piece = min(F_piece, BALL_R/self.dt) # prevent blow up
                    
                    f_x = F_piece*dx
                    f_y = F_piece*dy
                    
                    F_r[i] += f_x*np.cos(self.agents[i].th) + f_y*np.sin(self.agents[i].th)
                    F_r[j] += -f_x*np.cos(self.agents[j].th) + -f_y*np.sin(self.agents[j].th)
                    
                    r_cm_a = ((self.agents[i].balls[a][0] - self.agents[i].cm[0])**2 + (self.agents[i].balls[a][1] - self.agents[i].cm[1])**2)**0.5
                    r_cm_b = ((self.agents[j].balls[b][0] - self.agents[j].cm[0])**2 + (self.agents[j].balls[b][1] - self.agents[j].cm[1])**2)**0.5

                    
                    Tau[i] += (f_y*np.cos(self.agents[i].th) + -f_x*np.sin(self.agents[i].th))*r_cm_a
                    Tau[j] += (-f_y*np.cos(self.agents[j].th) + f_x*np.sin(self.agents[j].th))*r_cm_b
                    
                if r2 < (2*L)**2:
                    self.agents[i].has_neighbours = True
                    self.agents[j].has_neighbours = True

                else:
                    self.agents[i].has_neighbours = False
                    self.agents[j].has_neighbours = False  
        
    
    def run(self): 
        # setup animation
        if self.animation == 'on':
            ax = set_canvas(self.frame_lim)
            print 'on'

        details = {'r': BALL_R, 'frame_lim': FRAME_LIM, 'run_time': self.run_time, 'N': self.N, 'L': NUM_BALLS}
        cm_x_traj = []
        cm_y_traj = []
        th_traj = []
            
            
        # main loop
        while self.time < self.run_time:
        
            #Initialize F and tau
            F_r = np.zeros(self.N) #radial force
            Tau = np.zeros(self.N) #torque
            
            # loop over all rod-rod interactions to compute F and tau
            for i in xrange(self.N):
                for j in xrange(i+1, self.N):
                    
                    
                    if ((self.agents[i].cm[0] - self.agents[j].cm[0])**2 - (self.agents[i].cm[0] - self.agents[j].cm[0])**2)**0.5 > BALL_R*2*DIVISION_LENGTH:
                        # break back to start of for loop, skipping the ball-ball computations that follow
                        continue
                    
                    self.compute_forces(i, j, F_r, Tau)
                        
                        
            #update all agents
            for i in xrange(self.N):
                self.agents[i].step(F_r[i], Tau[i])
                
                cm_x_traj.append(self.agents[i].cm[0])
                cm_y_traj.append(self.agents[i].cm[1])
                th_traj.append(self.agents[i].th)
                
                # growth
                if self.agents[i].t > self.grow_time and self.agents[i].has_neighbours == False:
                    self.agents[i].N += 1
                    self.agents[i].balls = self.agents[i].compute_rod()
                    
                    self.agents[i].t = 0
                    
                    #division
                    if self.agents[i].N == self.division_length:

                        new_cm1 = self.agents[i].balls[int(np.floor(NUM_BALLS/2.))]
                        new_cm2 = self.agents[i].balls[self.division_length - int(np.floor(NUM_BALLS/2.)) -1]

                        self.agents[i].cm = new_cm1
                        self.agents[i].N = NUM_BALLS
                        
                        new_bacteria = Agent(cm = new_cm2, v0 = 0, ball_r = BALL_R, theta = self.agents[i].th, num_balls = NUM_BALLS, frame_lim = FRAME_LIM)
                        self.agents.append(new_bacteria)
                        
                        self.agents[i].balls = self.agents[i].compute_rod()
                        self.agents[-1].balls = self.agents[-1].compute_rod()
                    
                        self.N += 1
            
            if self.animation == 'on':      
                if self.time % self.show_every == 0:
                    draw_(ax, self.agents)
            
            self.time += 1
        
        
        return [cm_x_traj, cm_y_traj, th_traj], details

    