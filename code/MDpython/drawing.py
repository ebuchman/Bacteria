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
    
    
def compute_rod(cm_x, cm_y, th, N, r, frame_lim):
    balls = []
    for i in xrange(N):
        x = (cm_x - (N - abs(N%2 - 1) - 2*i)*r*np.cos(th))%frame_lim
        y = (cm_y - (N - abs(N%2 - 1) - 2*i)*r*np.sin(th))%frame_lim
        
        balls.append([x, y])
    return balls
    
def draw__(ax, balls, r):
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
    
def animation():
    import cPickle as pickle
    
    
    f = open('trajectories.pkl', 'rb')
    [cm_x_traj, cm_y_traj, th_traj], details = pickle.load(f)
    f.close()
    
    r = details['r']
    run_time = details['run_time']
    frame_lim = details['frame_lim']
    N = details['N']
    L = details['L']
    
    
    ax = set_canvas(frame_lim)

    for i in xrange(run_time):
        balls = []
        for j in xrange(N):
            balls += compute_rod(cm_x_traj[i*N + j], cm_y_traj[i*N + j], th_traj[i*N + j], L, r, frame_lim)
        
        if i%10 == 0:
            draw__(ax, balls, r)

if __name__ == "__main__":
    animation()