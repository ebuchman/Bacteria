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
  //*idum = (long) 271; // consistent
  *idum = -(long) time(NULL); // random 

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
	  struct Forces *fint, double dt, int t)
{
    
  double dvx, dvy, domega, dx, dy, dth;
  double net_f, net_th;
  double fx, fy, tau;
  double rod_end_x, rod_end_y;
  double dxp, dyp, R, L;
  double dU, dxs;
  struct Pillus * pil;
  
  struct pilForces pil_forces; // holds pillae forces summed over pillae (single agent)

  /* ROB: 
     
  Extend Pilus subroutine - should this be separate from the subroutine that 
  moves the bacteria (i.e. stepper)?

  */

  
  int j;

  // extend pilli
  for (j = 0; j < agents[i].Npil; j++)
  {
    if (agents[i].pillae[j].L <= 0)
    {
      if (ran1(idum) < p.PROB_EXTEND)
      {
        printf("extending pillus!\n");
        extend_pillus(agents[i].pillae, j, agents[i], idum, p);
      }
    }
  }

  compute_pilli_forces(&pil_forces, agents, i, p);
  
  net_th = compute_new_angle(pil_forces.Fx, pil_forces.Fy);

  /* Net force */
  fx = fint->Fx[i] + pil_forces.Fx;
  fy = fint->Fy[i] + pil_forces.Fy;
      
  //FRICTION    
      
  // static: if net force is greater than friction, otherwise, forces are 0
  if (agents[i].vx == 0 && agents[i].vy == 0)
  {
      net_f = sqrt(fx*fx + fy*fy) - p.STATIC_FRICTION;
      
      if (net_f > 0)
      {
        fx = fx - p.STATIC_FRICTION*cos(net_th);
        fy = fy - p.STATIC_FRICTION*sin(net_th); 
      }
      else
      {
        fx = 0;
        fy = 0;
      }
  }
  //kinetic: if greater than friction, subtract friction.  otherwise, set to zero.  add velocity-dependent friction
  else
  {
      net_f = sqrt(fx*fx + fy*fy) - p.KINETIC_FRICTION;

      if (net_f > 0)
      {
        fx = fx - p.KINETIC_FRICTION*cos(net_th);
        fy = fy - p.KINETIC_FRICTION*sin(net_th); 
      }
      else
      {
        fx = 0;
        fy = 0;
      }
      
      // velocity dependent dissipation
      fx -= p.GAMMA*agents[i].vx;
      fy -= p.GAMMA*agents[i].vy;
  }

  printf("i, fx, fy: %d, %f, %f\n", i, fx, fy);
  
  /* Net torque */
  tau = fint->Tau[i] + pil_forces.Tau - p.GAMMA*agents[i].omega;
  
  
  /* 
     ROB: Ok, this is the velocity-verlet routine, which will not be 
     modified. This could be in its own subroutine. 
  */

  
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
      
//      printf("L and R: %f, %f\n", L, R);
      //printf("x, y: %f, %f\n", dxp, dyp);
      //printf("th1, th2: %f, %f\n", agents[i].th*180./M_PI, pil->th*180./M_PI);
      
      //printf("x: %f; dx dy: %f, %f\n", pil->x_ext, dxp, dyp);
      
      //if (fabs((agents[i].th - pil->th)*180./M_PI) > 30) printf("diff th: %f \n", (agents[i].th - pil->th)*180./M_PI);
      
      //printf("diff th: %f \n", (agents[i].th - pil->th)*180./M_PI);
      //printf("fx, fy, tau: %f, %f, %f\n", fx, fy, tau);

      /*
      printf("cmx, cmy: %f, %f; rodx, rody: %f, %f; pilx, pily: %f, %f\n", agents[i].cm_x, agents[i].cm_y, rod_end_x, rod_end_y, pil->x, pil->y);
      
      printf("fx: %f, fy: %f L: %f, R: %f\n", fx, fy, L, R);
      printf("th: %f, dx, dy: %f, %f; dxp, dyp: %f, %f\n\n", agents[i].th, dx, dy, dxp, dyp);
    */
      if (R > L)
      {
        pil->x_ext += R - L; //extend sping
      }
      else if (R < L)
      {
        pil->L += R - L; //retract pillus
        //printf("\n retract, R, L, dx, dy: %f, %f, %f, %f, %f\n", pil->L, R, L, dx, dy);
      }
      
      dx = fabs(dx);
      dy = fabs(dy);
      
      if (dx < 0.0001 && dy < 0.0001) //|| R==L)
      {
        printf("no motion\n");
        pil->x_ext += (pil->P / p.STATIC_FRICTION)*p.DT;
      }
      

      pil->th = compute_new_angle(dxp, dyp);
      
      
      int snapped;
      
      
      if (pil->x_ext > 0.5*pil->L0) 
      //if (t == 10)
      { 
        printf("snapped!, %d\n\n", t);
        pil->L = 0.0; // pillus snaps     
        snapped = 1;
      }
      
      if (pil->L <= 0.0)
      {
        if (snapped != 1) {printf("done\n"); }//exit(0);}
        pil->L = 0.0;
        pil->x_ext = 0.0;
        pil->F = 0.0;   
        //exit(0);
      }

    }
    
    
    

  }
  
  
  //pbc           // ?????? 3 AND 6???
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
