import numpy as np
import matplotlib.pyplot as plt

from pylab import rand

RUN_P = 0.85
RUN_R = 0.05

TUMBLE_P = 0.3
TUMBLE_R = 20


BACTERIA_W = 0.2
BACTERIA_H = 0.6

class Agent():
    def __init__(self, xy, stiffness, frame_lim ):
        self.w = BACTERIA_W
        self.h = BACTERIA_H
        
        self.xy = xy
        self.th = np.random.rand()*360
        
        self.run_p = RUN_P#self.stiffness_function(stiffness)
        self.tumble_p = TUMBLE_P#self.run_p / 3.

        print self.th, self.run_p


        self.TUMBLE = False
        
        self.frame_lim = frame_lim
        
    
    def stiffness_function(self, stiffness):
        #logistic
        p = 1./(1 + np.e**(-stiffness))
        
        return p
    
    def step(self):
        # if not tumbling, pick random number.  If less than RUN_P, move RUN_R in direction of orientation.  else, start tumbling.
        # if tumbling, pick random number.  If greater than TUMBLE_P, rotate by TUMBLE_R.  else, stop tumbling.
        # matplotlib has (0,0) in the upper left - adapt trig accordingly...
                
        if not self.TUMBLE:
            p = rand()
            
            if p < self.run_p:
                self.xy[0] = (self.xy[0] + RUN_R*np.sin(self.th*np.pi/180.))%self.frame_lim
                self.xy[1] = (self.xy[1] - RUN_R*np.cos(self.th*np.pi/180.))%self.frame_lim
            else:
                self.TUMBLE = True
        
        if self.TUMBLE:
            p = rand()
            
            if p > self.tumble_p:
                q = rand()
                if q > 0.5:
                    self.th = (self.th + TUMBLE_R)%360
                else:
                    self.th = (self.th - TUMBLE_R)%360
            else:
                self.TUMBLE = False
        

