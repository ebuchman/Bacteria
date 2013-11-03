#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"


void step(struct Parameters p, long *idum, struct Box * grid, struct Agent *agents, double dt, int t);
void evolution(struct Parameters p, long *idum, struct Agent *agents, char * path, struct Box *grid);
void verlet(double fx, double fy, double tau, struct Agent * agents, int i, double dt);

/*****************************************************************************/
/*
  Stepper function 
*/

void step(struct Parameters p, long *idum, struct Box *grid, struct Agent *agents, double dt, int t)
{
  int i;
  double fx, fy, tau;  
  
  for (i = 0; i <p.NUM_BACTERIA; i++)
  {
    /* Net force */
    fx = agents[i].iFx + agents[i].pFx + agents[i].fFx;
    fy = agents[i].iFy + agents[i].pFy + agents[i].fFy;
      
    /* Net torque.  add an fTau (and have the gamma term in there)... */
    tau = agents[i].iTau + agents[i].pTau - p.GAMMA*agents[i].omega;

    //VERLET
    verlet(fx, fy, tau, agents, i, dt);

    //pbc
    pbc_position(agents, i, p);

    compute_rod(p, agents, i);

    agents[i].last_Fx = fx;
    agents[i].last_Fy = fy;
    agents[i].last_tau = tau;
    
    if (p.GRID == 1) update_grid_position(p, i, agents, grid);

  }
}


/*****************************************************************************/
/*
  This subroutine performs the time evolution.
*/

void evolution(struct Parameters p, long *idum, struct Agent *agents, char * path, struct Box *grid)
{
  int count = 0, t = 0;
  
  output_clean(p, path);

  while (t <= p.RUN_TIME)
  { /* Protocol:
            compute total force: interaction, pilli, friction
            evolve positions (verlet step)
            update pilli (includes re extension) */
      
      if (p.GRID == 1)
        compute_forces_grid(p, agents, grid);
      else
        compute_forces(p, agents);
      
      compute_pilli_forces(agents, p);
      
      friction(agents, p);

      step(p, idum, grid, agents, p.DT, t);
      
      update_pilli(agents, p, t, idum, grid);
      
      /* Periodically store results in files */

      if (t%p.SKIP == 0) 
      {
          data_out(p, agents, path);
          multiple_out(p, agents, count, path);
          count++;
      }
      
      if (t%(p.RUN_TIME/100) == 0) printf("# run time = %d\n", t);
      
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
