import matplotlib
matplotlib.use('macosx')

import numpy as np

import matplotlib.pyplot as plt
from pylab import figure, show
from matplotlib.patches import Ellipse

# this script animates trajectories of agents generated in C.  


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


def compute_rod(cm_x, cm_y, th, L, r, frame_lim):
    balls = []
    for i in xrange(L):
        x = (cm_x - (L - abs(L%2 - 1) - 2*i)*r*np.cos(th))%frame_lim
        y = (cm_y - (L - abs(L%2 - 1) - 2*i)*r*np.sin(th))%frame_lim
        
        balls.append([x, y])
    return balls

def open_file(filename, mode = None, sep = None):
    f = open(filename, 'r')
    if mode == 'np':
        data = np.fromfile(f, sep=sep)
    else:
        data = f.read()
    f.close()
    return data
    
def animation():
    cm_x_traj = open_file('cm_x_traj.txt', mode='np', sep = ',')
    cm_y_traj = open_file('cm_y_traj.txt', mode='np', sep =',')
    th_traj = open_file('th_traj.txt', mode='np', sep =',')
    
    details = open_file('details.txt', sep=None)
    
    #import ast
    #details = ast.literal_eval(details)
    details = eval(details)
        
    frame_lim = np.cast['int'](details['frame_lim'])
    r = np.cast['double'](details['r'])
    run_time = np.cast['int'](details['run_time'])
    N = np.cast['int'](details['N'])
    L = np.cast['int'](details['L'])
    
    ax = set_canvas(frame_lim)
    
    for i in xrange(run_time):
        balls = []
        for j in xrange(N):
            balls += compute_rod(cm_x_traj[i*N + j], cm_y_traj[i*N + j], th_traj[i*N + j], L, r, frame_lim)
        
        if i%10 == 0:
            draw_(ax, balls, r)
    

animation()
    
    
    
    