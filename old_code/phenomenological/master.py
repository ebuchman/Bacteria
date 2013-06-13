import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, show
from matplotlib.patches import Ellipse

from run_and_tumble import Agent

FRAME_LIM = 5

BACTERIA_W = 0.2
BACTERIA_H = 0.6

class Flock():
    def __init__(self, num_agents, stiffness, run_time, show_every):
        self.N = num_agents
        self.run_time = run_time
        self.show_every = show_every
                
        self.frame_lim = FRAME_LIM
        
        self.time = 0
        
        self.agents = []
        for i in xrange(self.N):
            xy = np.random.rand(2)*FRAME_LIM
            new_agent = Agent(xy, stiffness, FRAME_LIM)
            self.agents.append(new_agent)
            
    def run(self):
        # display stuff
        plt.ion()
        fig = figure()
        ax = fig.add_subplot(111, aspect='equal')
        
        ax.set_xlim(0, self.frame_lim)
        ax.set_ylim(0, self.frame_lim)

        # main loop
        while self.time < self.run_time:
            
            for a in self.agents:
                a.step()
                
            for i in xrange(self.N):
                for j in xrange(i+1, self.N):
                    delta_th = (self.agents[i].th - self.agents[j].th)
                    delta_r = ((self.agents[i].xy[0] - self.agents[j].xy[0])**2 + (self.agents[i].xy[1] - self.agents[j].xy[1])**2)**0.5
                    F = 10*np.e**(-delta_r*2/BACTERIA_H)*np.sin(delta_th*np.pi/180.)          #*np.sign(np.sin(delta_th*np.pi/180.))
                    self.agents[i].th += F
                    self.agents[j].th -= F
                    
                    self.agents[i].xy[0] = (self.agents[i].xy[0] + -0.01*(max(1./(delta_r/BACTERIA_H)**2, 0))*np.sin(self.agents[i].th*np.pi/180.))%self.frame_lim
                    self.agents[i].xy[1] = (self.agents[i].xy[1] - -0.01*(max(1./(delta_r/BACTERIA_H)**2, 0))*np.cos(self.agents[i].th*np.pi/180.))%self.frame_lim
                    
                    self.agents[j].xy[0] = (self.agents[j].xy[0] + -0.01*(max(1./(delta_r/BACTERIA_H)**2, 0))*np.sin(self.agents[j].th*np.pi/180.))%self.frame_lim
                    self.agents[j].xy[1] = (self.agents[j].xy[1] - -0.01*(max(1./(delta_r/BACTERIA_H)**2, 0))*np.cos(self.agents[j].th*np.pi/180.))%self.frame_lim
                    
            if self.time % self.show_every == 0:
                ax.clear()
                
                ells = []
                for a in self.agents:
                    ells.append(Ellipse(xy=a.xy, width=a.w, height=a.h, angle=a.th))

                for e in ells:
                    ax.add_artist(e)
                    e.set_clip_box(ax.bbox)
                    
                    #e.set_alpha(rand())
                    #e.set_facecolor(rand(3))
                
                plt.draw()

            self.time += 1

culture = Flock(num_agents = 20, stiffness = .3, run_time = 1000, show_every = 1)
culture.run()

plt.ioff()
    