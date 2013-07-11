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

  struct Forces forces; // holds MD forces in array over all agents

  struct Agent *agents;

  struct Parameters p;

  char path[100];

  idum = (long *)malloc(sizeof(long));
  *idum = (long) 271; // consistent
  //*idum = -(long) time(NULL); // random 

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

      /* ROB: Is this line necessary? */
      // isnt it? //
      agents[i].pillae = (struct Pillus *)malloc(p.NPIL*sizeof(struct Pillus));
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
	  struct Forces *fint, double dt, int t)
{

  double fx, fy, tau;
  double vx, vy;
  
  struct pilForces pil_forces; // holds pillae forces summed over pillae (single agent)
  
  extend_pilli(p, idum, i, agents);
  
  compute_pilli_forces(&pil_forces, agents, i, p);
    
  /* Net force */
  fx = fint->Fx[i] + pil_forces.Fx;
  fy = fint->Fy[i] + pil_forces.Fy;

  /* Net torque, without friction, for now */
  tau = fint->Tau[i] + pil_forces.Tau - p.GAMMA*agents[i].omega;

      
  // FRICTION 
  friction(&fx, &fy, agents, i, p, pil_forces);
  
  //VERLET
  verlet(fx, fy, tau, agents, i, dt);
  
  /* Update the pilli */
  update_pilli(agents, i, p);
  
  
  //printf("xext, L0, L, before: %f, %f, %f\n", agents[0].pillae[0].x_ext, agents[0].pillae[0].L0, agents[0].pillae[0].L);
  printf("time: %d\n", t);
  
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
        step(p, idum, i, agents, forces, p.DT, t);
      }
      printf("\n\n");

      /* Periodically store results in files */

      if (t%p.SKIP == 0) 
      {
          data_out(p, agents, path);
          multiple_out(p, agents, count, path);
          count++;
      }

      if (t%(p.RUN_TIME/100) == 0) printf("# run time = %d\n", t);

      /* Sanity check */

      if (t == 10000) 
      {
        printf("anchor x and y: %f, %f\n", agents[0].pillae[0].x, agents[0].pillae[0].y);
        exit(0);
      }
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
