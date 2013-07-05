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
  *idum = (long) 2718; // consistent
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
    
  double dvx, dvy, domega, dx, dy, dth;

  double fx, fy, tau;
  double rod_end_x, rod_end_y;
  double dxp, dyp, R, L;
  double dU, dxs;
  struct Pillus * pil;

  
    
  int j;


  for (j = 0; j < agents[i].Npil; j++)
  {
    if (agents[i].pillae[j].L <= 0)
    {
      if (1) //(ran1(idum) < p.PROB_EXTEND)
      {
        printf("extending pillus!\n");
        extend_pillus(agents[i].pillae, j, agents[i], idum, p);
      }
    }
  }

  compute_pilli_forces(fint, agents, i, p);

  /* Net force */
  
  fx = fint->Fx[i] - p.GAMMA*agents[i].vx;
  fy = fint->Fy[i] - p.GAMMA*agents[i].vy;
  
  /* Net torque */
  
  tau = fint->Tau[i] - p.GAMMA*agents[i].omega;
  
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
  agents[i].th = agents[i].th + dth;
  
  /* Update the pilli */
  
  rod_end_x = agents[i].cm_x + p.BALL_R*p.BACTERIA_LENGTH*cos(agents[i].th);
  rod_end_y = agents[i].cm_y + p.BALL_R*p.BACTERIA_LENGTH*sin(agents[i].th);
  
  rod_end_x = fmod(rod_end_x, p.SCREEN_W);
  rod_end_y = fmod(rod_end_y, p.SCREEN_W);
  
  for (j=0; j < agents[i].Npil; j++)
  {
    
    pil = &agents[i].pillae[j];
    
    if (pil->L > 0)
    {
      dxp = min_sep(p, pil->x, rod_end_x);
      dyp = min_sep(p, pil->y, rod_end_y);
      R = sqrt(dxp*dxp + dyp*dyp);
      L = pil->L;
      /* If the new position is further from the anchor point, we extend the spring;'
            If it is closer, we retract the pillus.  */
      
      
      //printf("tau: %f\n\n", tau);
      
    
      printf("cmx, cmy: %f, %f; rodx, rody: %f, %f; pilx, pily: %f, %f\n", agents[i].cm_x, agents[i].cm_y, rod_end_x, rod_end_y, pil->x, pil->y);
      
      printf("fx: %f, fy: %f L: %f, R: %f\n", fx, fy, L, R);
      printf("th: %f, dx, dy: %f, %f; dxp, dyp: %f, %f\n\n", agents[i].th, dx, dy, dxp, dyp);
    
      if (R > L)
      {
        pil->x_ext += R - L; //extend sping
      }
      else if (R < L)
      {
        pil->L += R - L; //retract pillus
      }
      
      if ((dx == 0 && dy == 0) || (R == L))
      {
        pil->x_ext += (pil->P / p.F_FRICTION)*p.DT;
      }
      

      pil->th = compute_new_angle(dxp, dyp, agents[i].th);
      
      if (pil->x_ext > 0.5*pil->L0) 
      { 
        printf("snapped!, %d\n\n", t);
        pil->L = 0; // pillus snaps     
      }
      if (pil->L == 0.0)
      {
        
        pil->x_ext = 0;
        pil->F = 0;   
        //exit(0);
      }

    }
    
    
    

  }
  
  
  //pbc           // ?????? 3 AND 6???
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
      
      /* Evolve positions */
      for (i = 0; i < p.NUM_BACTERIA; i++)
      {
        step(p, idum, i, agents, forces, p.DT, t);
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

      if (t == 3000) exit(0);

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
