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

  /* ROB: SHOULD THIS BE p = load_params(&p) ? */

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

      /* ROB: Is this line necessary? */

      agents[i].pillae = (struct Pillus *)malloc(p.NPIL*sizeof(struct Pillus));
    }    

  /* Initialize colony variables */

  if (p.UNIFORM == 0)
    make_colony(p, agents, idum); 
  else
    make_colony_uniform(p, agents);
      
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

  double fx, fy, tau;
    
  int j;

  /* ROB: 
     
  Extend Pilus subroutine - should this be separate from the subroutine that 
  moves the bacteria (i.e. stepper)?

  */

  for (j = 0; j < agents[i].Npil; j++)
  {
    if (agents[i].pillae[j].L <= 0)
    {
      if (ran1(idum) < p.PROB_EXTEND)
      {
        extend_pillus(agents[i].pillae, j, agents[i], idum);
      }
    }
  }
  
  compute_pilli_forces(fint, agents, i, p);

  /* Net force */

  fx = fint->Fx[i] - p.GAMMA*agents[i].vx;
  fy = fint->Fy[i] - p.GAMMA*agents[i].vy;
  
  /* Net torque */

  /* 
     ROB: Ok, this is the velocity-verlet routine, which will not be 
     modified. This could be in its own subroutine. 
  */

  tau = fint->Tau[i] - p.GAMMA*agents[i].omega;
  
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
  agents[i].th = agents[i].th + dth;
  

  /* 
     ROB: Are you rotating the pili too? shouldn't these remain fixed in 
     orientation while the bacteria rotate?
  */
  
  for (j=0; j < agents[i].Npil; j++)
  {
    agents[i].pillae[j].th -= dth;
  
  }
  
  //pbc

  /* ROB: You could write a subroutine to implement pbc (and use this 
     all over the place)
  */

  agents[i].cm_x = fabs(fmod(agents[i].cm_x + 3.0*p.SCREEN_W, p.SCREEN_W));
  agents[i].cm_y = fabs(fmod(agents[i].cm_y + 3.0*p.SCREEN_W, p.SCREEN_W));

  agents[i].th = fabs(fmod(agents[i].th + 6.0*M_PI, 2.0*M_PI));

  compute_rod(p, agents[i].balls, agents[i].cm_x, agents[i].cm_y, 
	      agents[i].ball_r, agents[i].th, agents[i].N);

  agents[i].last_Fx = fx;
  agents[i].last_Fy = fy;

  agents[i].last_tau = tau;
  /*
  if (agents[i].vx * cos(agents[i].th) < 0)
  {
    //reverse direction if there was a bounce back
    agents[i].th -= M_PI;
    for (j=0; j<agents[i].Npil; i++)
    {
      agents[i].pillae[j].L = 0;
    }
  }*/
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
      
      /* 
	 Evolve positions 

	 ROB: Are you evolving the bacteria postions sequentially and using 
	 the updated ith position to compute the motion of the (i+1)th bacteria
	 OK, not if all your bacteria interactions are in the compute_forces.
      */

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
