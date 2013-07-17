#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"


/*****************************************************************************/
/*
  Stepper function 
*/

void step(struct Parameters p, long *idum, int i, struct Agent *agents, double dt, int t)
{

  double fx, fy, tau;
  double vx, vy;
  
  extend_pilli(p, idum, i, agents);
  
  compute_pilli_forces(agents, i, p);
  printf ("pfx, pfy, ptau: %f, %f, %f\n", agents[i].pFx, agents[i].pFy, agents[i].pTau);
    
  /* Net force */
  fx = agents[i].iFx + agents[i].pFx;
  fy = agents[i].iFy + agents[i].pFy;

  /* Net torque, without friction, for now */
  tau = agents[i].iTau + agents[i].pTau - p.GAMMA*agents[i].omega;
      
  // FRICTION 
  friction(&fx, &fy, agents, i, p);
  
  //VERLET
  verlet(fx, fy, tau, agents, i, dt);
  
  /* Update the pilli */
  update_pilli(agents, i, p);
  
  
  //printf("xext, L0, L, before: %f, %f, %f\n", agents[0].pillae[0].x_ext, agents[0].pillae[0].L0, agents[0].pillae[0].L);
  //printf("time: %d\n", t);
  
  //pbc
  /* ROB: You could write a subroutine to implement pbc (and use this 
     all over the place)
  */
  agents[i].cm_x = fabs(fmod(agents[i].cm_x + 3.0*p.SCREEN_W, p.SCREEN_W));
  agents[i].cm_y = fabs(fmod(agents[i].cm_y + 3.0*p.SCREEN_W, p.SCREEN_W));
  agents[i].th = fabs(fmod(agents[i].th + 6.0*M_PI, 2.0*M_PI));


  compute_rod(p, agents, i); //agents[i].balls, agents[i].cm_x, agents[i].cm_y, agents[i].ball_r, agents[i].th, agents[i].N);

  agents[i].last_Fx = fx;
  agents[i].last_Fy = fy;
  agents[i].last_tau = tau;
    
}


/*****************************************************************************/
/*
  This subroutine performs the time evolution.
*/

void evolution(struct Parameters p, long *idum, struct Agent *agents, char * path, struct Box *grid)
{
  int i, count, t = 0;

  output_clean(p, path);

  count = 0;

  while (t <= p.RUN_TIME)
    {
      /* Compute forces */
      
      if (p.GRID == 1)
        compute_forces_grid(p, agents, grid, p.DT);
      else
        compute_forces(p, agents, p.DT);
//      printf("forces on %d: %f, %f\n", i, agents[0].iFx, agents[0].iFy);
      
      /* Evolve positions */

      for (i = 0; i < p.NUM_BACTERIA; i++)
      {
        step(p, idum, i, agents, p.DT, t);
        if (p.GRID == 1) update_grid_position(p, i, agents, grid);
        
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
