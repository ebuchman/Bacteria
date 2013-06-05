import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, show
from matplotlib.patches import Ellipse

def set_canvas(frame_lim):
    plt.ion()
    fig = figure()
    ax = fig.add_subplot(111, aspect='equal')
    
    ax.set_xlim(0, frame_lim)
    ax.set_ylim(0, frame_lim)
    
    return ax


def draw_(ax, agents):
    ax.clear()
    
    balls = []
    for a in agents:
        for ball in a.balls:
            balls.append(Ellipse(xy=ball, width=a.r*2, height=a.r*2, angle=0))

    for b in balls:
        ax.add_artist(b)
        b.set_clip_box(ax.bbox)
        
        #e.set_alpha(rand())
        #e.set_facecolor(rand(3))
    
    plt.draw()