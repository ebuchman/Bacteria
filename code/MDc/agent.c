#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "agent.h"

int RUN_TIME = 1000;

int NUM_BACTERIA = 10;
int BACTERIA_LENGTH = 4;
double SCREEN_W = 10;
double BALL_R = 0.2;
double GAMMA = 1.;
double E = 1.;

double DT = 0.001;


struct Forces * compute_forces(struct Agent ** agents, double dt)
{
    double L = BALL_R*2;
    struct Forces * forces = malloc(sizeof(struct Forces));

    forces->F_r = malloc(NUM_BACTERIA*sizeof(double));
    forces->Tau = malloc(NUM_BACTERIA*sizeof(double));

    int i, j, a, b; // iterators
    double dx, dy, r2;
    double F_piece, f_x, f_y;
    double r_cm_a, r_cm_b;
        
        
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
                   
                    if (r2 < pow(pow(2, 1./6)*L, 2))
                    {
                        F_piece = (48*E)*(pow(L, 12)*pow(r2, -6) - 0.5*pow(L, 6)*pow(r2, -3))/r2;

                        if (F_piece > BALL_R/dt)
                            F_piece = BALL_R/dt;
                        
                        f_x = F_piece*dx;
                        f_y = F_piece*dy;
                        
                        forces->F_r[i] += f_x*cos(agents[i]->th) + f_y*sin(agents[i]->th);
                        forces->F_r[j] += -f_x*cos(agents[i]->th) + -f_y*sin(agents[i]->th);
                                                
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

struct Agent * make_agent(int N, double x, double y, double theta)
{
    struct Agent * ag = malloc(sizeof(struct Agent));

    ag->N = N;

    ag->cm_x = x;
    ag->cm_y = y;
    ag->th = theta;
    ag->omega = 0;
    
    ag->F_self = 1.;
    
    ag->ball_r = BALL_R; 
    
    ag->last_F = 0;
    ag->last_tau = 0;
    
    ag->t = 0;
    
    return ag;
}

struct Agent ** make_colony()
{
    struct Agent ** agents = malloc(NUM_BACTERIA); //array of bacteria 
    struct Agent * ag = malloc(sizeof(struct Agent));
    srand(time(NULL));

    double x, y, th;
    int i; //iterator
    for (i=0; i < NUM_BACTERIA; i++)
    {
        x = SCREEN_W*rand()/(double)RAND_MAX;
        y = SCREEN_W*rand()/(double)RAND_MAX;
        th = 2*3.14*rand()/(double)RAND_MAX;
        ag = make_agent(BACTERIA_LENGTH, x, y, th);
        agents[i] = ag;
    }
    return agents;
}


double * compute_rod(double cm_x, double cm_y, double r, double th, int N)
{
    double * balls = malloc(N*2*sizeof(double));
    
    int i;
    double x, y;
    for (i=0; i < N; i++)
    {
        x = fmod((cm_x - (N - abs(N%2 - 1) - 2*i)*r*cos(th)), SCREEN_W);
        y = fmod((cm_y - (N - abs(N%2 - 1) - 2*i)*r*sin(th)), SCREEN_W);
        balls[i*2] = x;
        balls[i*2+1] = y;
    }

    return balls;
}

void step(struct Agent * agent, double F, double tau, double dt)
{
    
    double dv, domega, dx, dy, dth;
    
    F = F + agent->F_self - GAMMA*agent->v;
    tau = tau + (rand()/(double)RAND_MAX - 0.5)*2*2*3.14 - GAMMA*agent->omega;

    dv = 0.5*(F + agent->last_F)*dt;
    domega = 0.5*(tau + agent->last_tau)*dt;
    
    agent->v += dv;
    agent->omega += domega;
    
    dx = agent->v*cos(agent->th)*dt + 0.5*F*cos(agent->th)*pow(dt, 2);
    dy = agent->v*sin(agent->th)*dt + 0.5*F*sin(agent->th)*pow(dt, 2);

    dth = agent->omega*dt + 0.5*tau*pow(dt, 2);
    
    
    agent->cm_x = fmod(agent->cm_x + dx, SCREEN_W);
    agent->cm_y = fmod(agent->cm_y + dy, SCREEN_W);
    agent->th = fmod(agent->th + dth, 2*3.14);
    
    agent->balls = compute_rod(agent->cm_x, agent->cm_y, agent->ball_r, agent->th, agent->N);
    
    agent->last_F = F;
    agent->last_tau = tau;
    
    
}

void compute_rod_all(struct Agent ** agents)
{
    int i;
    for (i = 0; i < NUM_BACTERIA; i++)
    {
        agents[i]->balls = compute_rod(agents[i]->cm_x, agents[i]->cm_y, agents[i]->ball_r, agents[i]->th, agents[i]->N); 
    
    }
}

int main()
{
    int i;
    int time = 0;
    double dt = DT;
    struct Forces * forces = malloc(sizeof(struct Forces));

    double * cm_x_traj = malloc(RUN_TIME*NUM_BACTERIA*sizeof(double));
    double * cm_y_traj = malloc(RUN_TIME*NUM_BACTERIA*sizeof(double));
    double * th_traj = malloc(RUN_TIME*NUM_BACTERIA*sizeof(double));
    

    struct Agent ** agents = make_colony();
    i=0;
    
    compute_rod_all(agents);
    
    while (time < RUN_TIME)
    {
        forces = compute_forces(agents, dt);
                    
        for (i = 0; i < NUM_BACTERIA; i++)
        {
            step(agents[i], forces->F_r[i], forces->Tau[i], dt);
            compute_rod_all(agents);
            
            cm_x_traj[time*NUM_BACTERIA + i] = agents[i]->cm_x;
            cm_y_traj[time*NUM_BACTERIA + i] = agents[i]->cm_y;
            th_traj[time*NUM_BACTERIA + i] = agents[i]->th;

        }
        time++;
    }
    
    
    free(agents);
    free(forces);
    
    printf("done computations.  saving trajectories...\n");
    
    FILE *fx = fopen("cm_x_traj.txt", "w");
    FILE *fy = fopen("cm_y_traj.txt", "w");
    FILE *fth = fopen("th_traj.txt", "w");

    for (i = 0; i < RUN_TIME*NUM_BACTERIA; i++)
    {
        fprintf(fx, "%f,",cm_x_traj[i]);
        fprintf(fy, "%f,",cm_y_traj[i]);
        fprintf(fth, "%f,",th_traj[i]);

    }
    fclose(fx);
    fclose(fy);
    fclose(fth);
    
    return 0;


} 