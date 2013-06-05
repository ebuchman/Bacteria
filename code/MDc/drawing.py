import matplotlib
matplotlib.use('macosx')

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


def draw_(ax, balls, r):
    ax.clear()
    
    full_balls = []
    for ball in balls:
        full_balls.append(Ellipse(xy=ball, width=r*2, height=r*2, angle=0))

    for b in full_balls:
        ax.add_artist(b)
        b.set_clip_box(ax.bbox)
        
        #e.set_alpha(rand())
        #e.set_facecolor(rand(3))
    
    plt.draw()


def compute_rod(cm_x, cm_y, th, N, r, frame_lim):
    balls = []
    for i in xrange(N):
        x = (cm_x - (N - abs(N%2 - 1) - 2*i)*r*np.cos(th))%frame_lim
        y = (cm_y - (N - abs(N%2 - 1) - 2*i)*r*np.sin(th))%frame_lim
        
        balls.append([x, y])
    return balls

def get_traj(filename):
    f = open(filename, 'r')
    traj = np.fromfile(f, sep=',')
    f.close()
    return traj

def animation():
    cm_x_traj = get_traj('cm_x_traj.txt')
    cm_y_traj = get_traj('cm_y_traj.txt')
    th_traj = get_traj('th_traj.txt')
    
    ax = set_canvas(10)
    
    for i in xrange(10000):
        balls = []
        for j in xrange(10):
            balls += compute_rod(cm_x_traj[i*10 + j], cm_y_traj[i*10 + j], th_traj[i*10 + j], 4, 0.2, 10)
        
        if i%50 == 0:
            draw_(ax, balls, 0.2)
    


animation()
    
    
    
    