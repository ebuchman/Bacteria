#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "agent.h"

/* Molecular Dynamics on self-propelled rods.  Main runs a single simulation, recording trajectories of all agents (cm, theta) in txt files, to be animated by drawing.py */

int RUN_TIME = 10000;

int NUM_BACTERIA = 10;
int BACTERIA_LENGTH = 4;
double SCREEN_W = 10;
double BALL_R = 0.2;
double GAMMA = 1.; //friction
double E = 1.; //energy scale
double PI = 3.141592653589793;


double DT = 0.01;

// return vector of forces and vector of torques for entire colony
struct Forces * compute_forces(struct Agent ** agents, double dt)
{
    double L = BALL_R*2;
    
    struct Forces * forces = malloc(sizeof(struct Forces)); //to return F_r and Tau

    forces->F_r = malloc(NUM_BACTERIA*sizeof(double));
    forces->F_t = malloc(NUM_BACTERIA*sizeof(double));
    forces->Tau = malloc(NUM_BACTERIA*sizeof(double));

    int i, j, a, b; // iterators
    double dx, dy, r2;
    double F_piece, f_x, f_y;
    double r_cm_a, r_cm_b; //cm's of balls
        
        
    for (i = 0; i < NUM_BACTERIA; i++)
    {
        for (j = i + 1; j < NUM_BACTERIA; j++)
        {
            for (a = 0; a < BACTERIA_LENGTH; a ++)
            {
                for (b = 0; b < BACTERIA_LENGTH; b++)
                {
                    dx = agents[i]->balls[a*2] - agents[j]->balls[b*2];
                    dy = agents[i]->balls[a*2 + 1] - agents[j]->balls[b*2 + 1];
                    r2 = pow(dx, 2) + pow(dy, 2);
                    
                    if (r2 < pow(pow(2, 1./6)*L, 2)) //if close enough ...
                    {
                        F_piece = (48*E)*(pow(L, 12)*pow(r2, -6) - 0.5*pow(L, 6)*pow(r2, -3))/r2;

                        if (F_piece > BALL_R/dt) // cap it for stability ...
                            F_piece = BALL_R/dt; 
                        
                        f_x = F_piece*dx;
                        f_y = F_piece*dy;
                        
                        forces->F_r[i] += f_x*cos(agents[i]->th) + f_y*sin(agents[i]->th);
                        forces->F_r[j] += -f_x*cos(agents[j]->th) + -f_y*sin(agents[j]->th);
                        
                        forces->F_t[i] += -f_x*sin(agents[i]->th) + f_y*cos(agents[i]->th);
                        forces->F_t[j] += f_x*sin(agents[j]->th) - f_y*cos(agents[j]->th);
                        
                                                
                        r_cm_a = pow(pow(agents[i]->balls[a*2] - agents[i]->cm_x, 2) + pow(agents[i]->balls[a*2+1] - agents[i]->cm_y, 2), 0.5);
                        r_cm_b = pow(pow(agents[j]->balls[b*2] - agents[j]->cm_x, 2) + pow(agents[j]->balls[b*2+1] - agents[j]->cm_y, 2), 0.5);
                        
                        forces->Tau[i] += (f_y*cos(agents[i]->th) + -f_x*sin(agents[i]->th))*r_cm_a;
                        forces->Tau[j] += (-f_y*cos(agents[j]->th) + f_x*sin(agents[j]->th))*r_cm_b;
                    }
                    
                }
            }
        }
    }
    
    return forces;
}

//extend an agent from his cm in direction of theta 
double * compute_rod(double cm_x, double cm_y, double r, double th, int N)
{
    double * balls = malloc(N*2*sizeof(double));
    
    int i;
    double x, y;
    double end = -(N-1); // number of radi from cm the rod extends
    
    for (i=0; i < N; i++)
    {
        //pbc
        
        x = fabs(fmod((cm_x + (end + 2*i)*r*cos(th)), SCREEN_W));
        y = fabs(fmod((cm_y + (end + 2*i)*r*sin(th)), SCREEN_W));
        balls[i*2] = x;
        balls[i*2+1] = y;
    }

    return balls;
}

//convenient
void compute_rod_all(struct Agent ** agents)
{
    double * balls = malloc(BACTERIA_LENGTH*2*sizeof(double));
    int i;
    for (i = 0; i < NUM_BACTERIA; i++)
    {
        balls = compute_rod(agents[i]->cm_x, agents[i]->cm_y, agents[i]->ball_r, agents[i]->th, agents[i]->N);
        agents[i]->balls =  balls;
    
    }
}


// return pointer to new agent
struct Agent * make_agent(int N, double x, double y, double theta)
{
    struct Agent * ag = malloc(sizeof(struct Agent));

    ag->N = N;

    ag->cm_x = x;
    ag->cm_y = y;
    ag->th = theta;
    ag->omega = 0;
    ag->v_r = 0;
    ag->v_t = 0;
    
    ag->F_self = 1.;
    
    ag->ball_r = BALL_R; 
    
    ag->last_F_r = 0;
    ag->last_F_t = 0;
    ag->last_tau = 0;
    
    ag->t = 0;
    
    return ag;
}

// return list of pointers to agents ... the pointer vs. array subtlety is bugging me a bit ...
struct Agent ** make_colony()
{
    struct Agent ** agents = malloc(NUM_BACTERIA); //array of bacteria 
    struct Agent * ag = malloc(sizeof(struct Agent)); //single bacteria
    

    //srand(time(NULL));

    double x, y, th;
    int i; //iterator
    for (i=0; i < NUM_BACTERIA; i++)
    {
        x = SCREEN_W*rand()/(double)RAND_MAX;
        y = SCREEN_W*rand()/(double)RAND_MAX;
        th = 2*PI*rand()/(double)RAND_MAX;
        
        // if I try to do these two lines in one line, eliminating the Agent * ag above, the thing breaks.  I'm curious to understand why... 
        ag = make_agent(BACTERIA_LENGTH, x, y, th);
        agents[i] = ag;
    }
    
    compute_rod_all(agents);    
    return agents;
}

//step
void step(struct Agent * agent, double F_r, double F_t, double tau, double dt)
{
    
    double dv_r, dv_t, domega, dx, dy, dth;
    double * balls = malloc(BACTERIA_LENGTH*2*sizeof(double));
    
    F_r = F_r + agent->F_self - GAMMA*agent->v_r;
    F_t = F_t - GAMMA*agent->v_t;
    tau = tau + (rand()/(double)RAND_MAX - 0.5)*2*2*PI - GAMMA*agent->omega;

    dv_r = 0.5*(F_r + agent->last_F_r)*dt;
    dv_t = 0.5*(F_t + agent->last_F_t)*dt;
    domega = 0.5*(tau + agent->last_tau)*dt;

    agent->v_r += dv_r;
    agent->v_t += dv_t;
    agent->omega += domega;
    
    dx = (agent->v_r*cos(agent->th) - agent->v_t*sin(agent->th))*dt + 0.5*(F_r*cos(agent->th) - F_t*sin(agent->th))*pow(dt, 2);
    dy = (agent->v_r*sin(agent->th) + agent->v_t*cos(agent->th))*dt + 0.5*(F_r*sin(agent->th) + F_t*cos(agent->th))*pow(dt, 2);

    dth = agent->omega*dt + 0.5*tau*pow(dt, 2);
    
    //pbc
    agent->cm_x = fabs(fmod(agent->cm_x + dx, SCREEN_W));
    agent->cm_y = fabs(fmod(agent->cm_y + dy, SCREEN_W));
    agent->th = fabs(fmod(agent->th + dth, 2*PI));
    
    balls =  compute_rod(agent->cm_x, agent->cm_y, agent->ball_r, agent->th, agent->N);
    agent->balls = balls; //compute_rod(agent->cm_x, agent->cm_y, agent->ball_r, agent->th, agent->N);
    
    agent->last_F_r = F_r;
    agent->last_F_t = F_t;
    agent->last_tau = tau;
    
    
}

// run a single simulation.  save trajectories of cm's and theta.
int main()
{

    int i;
    int time = 0;
    double dt = DT;
    struct Forces * forces = malloc(sizeof(struct Forces));



    /*trajectories to be written to txt, and animated from python.  I haven't investigated enough, but memory doesn't seem to like what happens when RUN_TIME*NUM_BACTERIA gets big.  Other approaches      to recording trajectory data? */
    double * cm_x_traj = malloc(RUN_TIME*NUM_BACTERIA*sizeof(double));
    double * cm_y_traj = malloc(RUN_TIME*NUM_BACTERIA*sizeof(double));
    double * th_traj = malloc(RUN_TIME*NUM_BACTERIA*sizeof(double));


    struct Agent ** agents = make_colony(); //initialize colony
            
    // run simulation
    while (time < RUN_TIME)
    {
        forces = compute_forces(agents, dt);
                    
        for (i = 0; i < NUM_BACTERIA; i++)
        {
            step(agents[i], forces->F_r[i], forces->F_t[i], forces->Tau[i], dt);
            cm_x_traj[time*NUM_BACTERIA + i] = agents[i]->cm_x;
            cm_y_traj[time*NUM_BACTERIA + i] = agents[i]->cm_y;
            th_traj[time*NUM_BACTERIA + i] = agents[i]->th;
        }
        
        if (time == 1)
        {
            printf("(%f, %f), %f\n\n", agents[0]->cm_x, agents[0]->cm_y, agents[0]->th);
            
            for(i = 0; i <BACTERIA_LENGTH; i++)
            {
                printf("(%f,%f)\t", agents[0]->balls[2*i], agents[0]->balls[2*i+1]);
            }
                printf("\n\n");

        }

        time++;

    }
    

    
    // legit?  I imagine there's a more ethical way to do this...
    free(agents);
    free(forces);
    
    printf("done computations.  saving trajectories...\n");
    
    //save data to text files for animation with python
    FILE *fd = fopen("details.txt", "w");
    FILE *fx = fopen("cm_x_traj.txt", "w");
    FILE *fy = fopen("cm_y_traj.txt", "w");
    FILE *fth = fopen("th_traj.txt", "w");

    // should I have opened the files at the beginning, and run this throughout?
    for (i = 0; i < RUN_TIME*NUM_BACTERIA; i++)
    {
        fprintf(fx, "%f,",cm_x_traj[i]);
        fprintf(fy, "%f,",cm_y_traj[i]);
        fprintf(fth, "%f,",th_traj[i]);

    }
    
    // print a python dictionary of simulation details
    fprintf(fd, "{'r': %f, 'frame_lim': %f, 'run_time': %d, 'N': %d, 'L': %d}", BALL_R, SCREEN_W, RUN_TIME, NUM_BACTERIA, BACTERIA_LENGTH);


    fclose(fd);
    fclose(fx);
    fclose(fy);
    fclose(fth);
    
    printf("toodles!\n");
    
    return 0;
} 