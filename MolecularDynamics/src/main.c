#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

int main()
{
  int i;

  long *idum;

  struct Forces forces;

  struct Agent *agents;

  struct Parameters p;

  char path[100];

  idum = (long *)malloc(sizeof(long));
  *idum = -(long) time(NULL);

  /* Load Parameters */
  p = load_params(p);
  
  /* Build Directory Structure */
  sprintf(path, "data/N%d_T%d/", p.NUM_BACTERIA, p.RUN_TIME);
  printf("here: %s\n", path);
  mk_dirs(path); // make directories for data

  /* Allocate memory */
  forces.Fx = (double *) malloc(p.NUM_BACTERIA*sizeof(double));
  forces.Fy = (double *) malloc(p.NUM_BACTERIA*sizeof(double));
  forces.Tau = (double *) malloc(p.NUM_BACTERIA*sizeof(double));
  
  /* Zero initial forces */
  for(i=0; i<p.NUM_BACTERIA; i++)
    {
      forces.Fx[i] = 0.0;
      forces.Fy[i] = 0.0;
      forces.Tau[i] = 0.0;
    }

  /* Initialize agents */
  agents = (struct Agent *)malloc(sizeof(struct Agent)*p.NUM_BACTERIA);

  for(i=0; i<p.NUM_BACTERIA; i++)
    {
      agents[i].balls = (double *)malloc(p.BACTERIA_LENGTH*2*sizeof(double));
    }    

  /* Initialize colony variables */
  make_colony(p, agents, idum);
    
  /*** RUN SIMULATION ***/
  evolution(p, idum, &forces, agents, path);

  printf("# toodles!\n");

  free(idum);
    
  free(agents);

  return 0;

}


/*****************************************************************************/
/*
  Stepper function 
*/

void step(struct Parameters p, long *idum, int i, struct Agent *agents, 
	  struct Forces *fint, double dt)
{
    
  double dvx, dvy, domega, dx, dy, dth;
  double rand_tau;
  double fx, fy, tau;
    
  int j;

  /* Net force */
  
  fx = fint->Fx[i] + p.F_SELF*cos(agents[i].th) - p.GAMMA*agents[i].vx;
  fy = fint->Fy[i] + p.F_SELF*sin(agents[i].th) - p.GAMMA*agents[i].vy;
  
  /* Net torque */
  
  rand_tau = 6*M_PI*ran1(idum) -3*M_PI);
  
  tau = fint->Tau[i] + rand_tau - p.GAMMA*agents[i].omega;

  /* Update the agent */
  
  dvx = 0.5*(fx + agents[i].last_Fx)*dt;
  dvy = 0.5*(fy + agents[i].last_Fy)*dt;
  domega = 0.5*(tau + agents[i].last_tau)*dt;
  
  agents[i].vx += dvx; 
  agents[i].vy += dvy;
  agents[i].omega += domega;

  dx = agents[i].vx*dt + 0.5*fx*dt*dt;
  dy = agents[i].vy*dt + 0.5*fy*dt*dt;
  dth = agents[i].omega*dt + 0.5*tau*dt*dt;
  
  agents[i].cm_x = agents[i].cm_x + dx;
  agents[i].cm_y = agents[i].cm_y + dy;
  
  // if agents are extended, adjust theta. Otherwise, theta is in direction of velocity
  if (p.BACTERIA_LENGTH > 1)
    agents[i].th = agents[i].th + dth;
  else
    agents[i].th = atan(agents[i].vy/agents[i].vx);
  
    
  //pbc           // ?????? 3 AND 6???
  agents[i].cm_x = fabs(fmod(agents[i].cm_x + 3.0*p.SCREEN_W, p.SCREEN_W));
  agents[i].cm_y = fabs(fmod(agents[i].cm_y + 3.0*p.SCREEN_W, p.SCREEN_W));

  agents[i].th = fabs(fmod(agents[i].th + 6.0*M_PI, 2.0*M_PI));

  compute_rod(p, agents[i].balls, agents[i].cm_x, agents[i].cm_y, 
	      agents[i].ball_r, agents[i].th, agents[i].N);

  agents[i].last_Fx = fx;
  agents[i].last_Fy = fy;

  agents[i].last_tau = tau;
  
}


/*****************************************************************************/
/*
  This subroutine performs the time evolution.
*/

void evolution(struct Parameters p, long *idum, struct Forces *forces, 
	       struct Agent *agents, char * path)
{
  int i, count, t = 0;

  output_clean(p, path);
  count = 0;



  while (t <= p.RUN_TIME)
    {
      /* Compute forces */
      compute_forces(p, forces, agents, p.DT);
      
      /* Evolve positions */
      for (i = 0; i < p.NUM_BACTERIA; i++)
      {
        step(p, idum, i, agents, forces, p.DT);
      }

      /* Periodically store results in files */

      if (t%p.SKIP == 0) 
      {
          data_out(p, agents, path);
          multiple_out(p, agents, count, path);
          count++;
      }

      if (t%(p.RUN_TIME/100) == 0) printf("# run time = %d\n", t);

      /* Sanity check */

      if (t == 1)
        {
	  printf("# (%f, %f), %f\n\n", agents[0].cm_x, agents[0].cm_y, 
		 agents[0].th);
            
	  for(i = 0; i <p.BACTERIA_LENGTH; i++)
            {
	      printf("# (%f,%f)\t", agents[0].balls[2*i], agents[0].balls[2*i+1]);
            }

	  printf("\n\n"); 
        }

      t++;

      
    }
}
