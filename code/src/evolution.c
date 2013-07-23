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

void step(struct Parameters p, long *idum, int i, struct Box *grid, struct Agent *agents, double dt, int t)
{

  double fx, fy, tau;
  double vx, vy;
  
  extend_pilli(p, idum, grid, agents, i);

  friction(agents, i, p);

  /* Net force */
  fx = agents[i].iFx + agents[i].pFx + agents[i].fFx;
  fy = agents[i].iFy + agents[i].pFy + agents[i].fFy;
    
  /* Net torque, without friction, for now */
  tau = agents[i].iTau + agents[i].pTau - p.GAMMA*agents[i].omega;

  //VERLET
  verlet(fx, fy, tau, agents, i, dt);
  
  /* Update the pilli */
  update_pilli(agents, i, p, t);
  
  //pbc
  pbc_position(agents, i, p);

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
        compute_forces_grid(p, agents, grid);
      else
        compute_forces(p, agents);

      compute_pilli_forces(agents, p);


//      printf("forces on %d: %f, %f\n", i, agents[0].iFx, agents[0].iFy);
      
      /* Evolve positions */

      for (i = 0; i < p.NUM_BACTERIA; i++)
      {
        //printf("time: %d\n", t);
        step(p, idum, i, grid, agents, p.DT, t);
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

void verlet(double fx, double fy, double tau, struct Agent * agents, int i, double dt)
{
    double dvx, dvy, domega, dx, dy, dth;

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
}
