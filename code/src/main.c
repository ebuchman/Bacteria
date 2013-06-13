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

  idum = (long *)malloc(sizeof(long));
  *idum = -(long) time(NULL);

  p = load_params(p);

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
      agents[i].pillae = (struct Pillus *)malloc(p.NPIL*sizeof(struct Pillus));
    }    

  /* Initialize colony variables */
  
  
  if (p.UNIFORM == 0) 
    make_colony(p, agents, idum); 
  else
    make_colony_uniform(p, agents);
    
  printf("# %d\t%lf\t%lf\n", agents[0].N, agents[0].cm_x, agents[0].cm_y);

  /*** RUN SIMULATION ***/
  evolution(p, idum, &forces, agents);

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

  double fx, fy, tau;
  
  int j;

  for (j = 0; j < agents[i].Npil; j++)
  {
    if (agents[i].pillae[j].L <= 0)
    {
      if (ran1(idum) > 0.8)  extend_pillus(agents[i].pillae, j, agents[i], idum);
    }
  }
  compute_pilli_forces(fint, agents, i, p);
    
  /* Net force */

  fx = fint->Fx[i] +  agents[i].F_self*cos(agents[i].th) - p.GAMMA*agents[i].vx;
  fy = fint->Fy[i] +  agents[i].F_self*sin(agents[i].th) - p.GAMMA*agents[i].vy;
  
  /* Net torque */

  tau = fint->Tau[i] + (ran1(idum) - 0.5)*2.0*2.0*M_PI - p.GAMMA*agents[i].omega;

  dvx = 0.5*(fx + agents[i].last_Fx)*dt;
  dvy = 0.5*(fy + agents[i].last_Fy)*dt;

  domega = 0.5*(tau + agents[i].last_tau)*dt;
  
  agents[i].vx += dvx; 
  agents[i].vy += dvy;
  
  if (agents[i].vx * cos(agents[i].th) < 0)
  {
    //reverse direction if there was a bounce back
    agents[i].th -= M_PI;
  }

  agents[i].omega += domega;

  dx = agents[i].vx*dt + 0.5*fx*dt*dt;
  dy = agents[i].vy*dt + 0.5*fy*dt*dt;

  dth = agents[i].omega*dt + 0.5*tau*dt*dt;

  agents[i].cm_x = agents[i].cm_x + dx;
  agents[i].cm_y = agents[i].cm_y + dy;
  agents[i].th = agents[i].th + dth;
    
  //pbc
  
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
	       struct Agent *agents)
{
  int i, count, t = 0;

  output_clean(p);
    
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
          data_out(p, agents);
          multiple_out(p, agents, count);

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
