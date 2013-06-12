#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "rob.h"

int main()
{
  int i;

  long *idum;

  struct Forces forces;

  struct Agent *agents;

  struct Parameters p;

  char dummy[100];
 
  FILE *fp;

  idum = (long *)malloc(sizeof(long));
  *idum = -(long) time(NULL);

  /* Read in parameters from input file */

  fp = fopen("input","r");

  fscanf(fp,"%s", dummy);

  fscanf(fp, "%s %d", dummy, &(p.RUN_TIME) );
  fscanf(fp, "%s %d", dummy, &(p.SKIP) );
  fscanf(fp, "%s %d", dummy, &(p.NUM_BACTERIA) );
  fscanf(fp, "%s %d", dummy, &(p.BACTERIA_LENGTH) );

  fscanf(fp, "%s %d", dummy, &(p.UNIFORM) );
  
  fscanf(fp, "%s %lf", dummy, &(p.SCREEN_W) );
  fscanf(fp, "%s %lf", dummy, &(p.BALL_R) );
  fscanf(fp, "%s %lf", dummy, &(p.GAMMA) );
  fscanf(fp, "%s %lf", dummy, &(p.E) );

  fscanf(fp, "%s %lf", dummy, &(p.DT) );

  fclose(fp);

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


void make_colony(struct Parameters p, struct Agent *agents, long *idum)
{ 
  double  x, y, th;
  int i;

  for (i=0; i < p.NUM_BACTERIA; i++)
    {
      
      x = p.SCREEN_W*ran1(idum);
      y = p.SCREEN_W*ran1(idum);

      th = 2.0*M_PI*ran1(idum);
      
      agents[i].N = p.BACTERIA_LENGTH;
  
      agents[i].cm_x = x;
      agents[i].cm_y = y;

      agents[i].th = th;
      
      agents[i].omega = 0;
      agents[i].vx = 0;
      agents[i].vy = 0;
      
      agents[i].F_self = 1.0;
      
      agents[i].ball_r = p.BALL_R; 
      
      agents[i].last_Fx = 0;
      agents[i].last_Fy = 0;

      agents[i].last_tau = 0;
      
      agents[i].t = 0;

      compute_rod(p, agents[i].balls, agents[i].cm_x, agents[i].cm_y, 
		  agents[i].ball_r, agents[i].th, agents[i].N);
    }  
}

/*****************************************************************************/


void make_colony_uniform(struct Parameters p, struct Agent *agents)
{ 
  double dx, dy, offx, offy, x, y, th;
  int i, slip;

  dx = 1.02*p.BACTERIA_LENGTH*p.BALL_R*2.0;
  dy = 1.1*p.BALL_R*2.0;
  
  offx = 0.001;
  offy = 0.001;

  x = offx;
  y = offy;

  slip = 0;

  for (i=0; i < p.NUM_BACTERIA; i++)
    { 
      if (x >= p.SCREEN_W + offx - dx) 
	{
	  y += dy;
	  
	  offx = 0.001;

	  if (slip == 0)
	    {
	      offx = dx/2.0 + 0.001;
	      
	      slip = 1;
	    }
	  else
	    {
	      slip = 0;
	    }

	  x = offx;
	  
	  if (y >= p.SCREEN_W + offy - dy) 
	    {
	      printf("# Error!\n");
	      y = fabs(fmod(y,p.SCREEN_W));
	    }
	}
      
      th = 0.0;

      /* This is useful to set the proper colony size */

      printf("%lf\t%lf\n", x, y);

      agents[i].N = p.BACTERIA_LENGTH;
  
      agents[i].cm_x = x;
      agents[i].cm_y = y;

      agents[i].th = th;
      
      agents[i].omega = 0;
      agents[i].vx = 0;
      agents[i].vy = 0;
      
      agents[i].F_self = 1.0;
      
      agents[i].ball_r = p.BALL_R; 
      
      agents[i].last_Fx = 0;
      agents[i].last_Fy = 0;

      agents[i].last_tau = 0;
      
      agents[i].t = 0;

      compute_rod(p, agents[i].balls, agents[i].cm_x, agents[i].cm_y, 
		  agents[i].ball_r, agents[i].th, agents[i].N);

      x += dx;

    }  
}

/*****************************************************************************/
/*
  extend an agent from his cm in direction of theta 
*/

void compute_rod(struct Parameters p, double *balls, double cm_x, double cm_y, 
		 double r, double th, int N)
{
  int i;
  double x, y;

  for (i=0; i<N; i++)
    {

      x = cm_x - (N - 1 - 2*i)*r*cos(th);
      y = cm_y - (N - 1 - 2*i)*r*sin(th);

      // pbc

      x = fabs(fmod(x + 3.0*p.SCREEN_W, p.SCREEN_W));
      y = fabs(fmod(y + 3.0*p.SCREEN_W, p.SCREEN_W));

      balls[i*2] = x;
      balls[i*2+1] = y;
    }
}

double min_sep(struct Parameters p, double a, double b)
{
  double ds;

  ds = a - b;

  if (ds > 0.5*p.SCREEN_W)
    {
      ds -= p.SCREEN_W;
    }
  else if (ds < -0.5*p.SCREEN_W)
    {
      ds += p.SCREEN_W;
    }

  return ds;
}
  

/*****************************************************************************/

// return vector of forces and vector of torques for entire colony

void compute_forces(struct Parameters p, struct Forces *forces, 
		    struct Agent *agents, double dt)
{
  double L = p.BALL_R*2;
    

  int i, j, a, b; // iterators
  
  double dx, dy, r2;
  double F_piece, f_x, f_y;
  double r_cm_a, r_cm_b; //cm displacement's of balls

        
  /* zero the forces */
    
  for (i = 0; i < p.NUM_BACTERIA; i++)
    {
      forces->Fx[i] = 0.0;
      forces->Fy[i] = 0.0;
      forces->Tau[i] = 0.0;
    }
  
  /* compute the new forces */

  for (i = 0; i < p.NUM_BACTERIA; i++)
    {
      for (j = i + 1; j < p.NUM_BACTERIA; j++)
        {

	  for (a = 0; a < p.BACTERIA_LENGTH; a ++)
            {
	      for (b = 0; b < p.BACTERIA_LENGTH; b++)
                {
		  
		  dx = min_sep(p, agents[i].balls[a*2], agents[j].balls[b*2]);
		  
		  dy = min_sep(p, agents[i].balls[a*2 + 1], agents[j].balls[b*2 + 1]);

		  r2 = dx*dx + dy*dy;
  
		  if (r2 < pow(pow(2, 1./6)*L, 2)) //if close enough ...
                    {
		      F_piece = (48*p.E)*(pow(L, 12)*pow(r2, -6) - 
					  0.5*pow(L, 6)*pow(r2, -3))/r2;

		    
		      if (F_piece > p.BALL_R/dt) // cap it for stability ...
		      	F_piece = p.BALL_R/dt; 

		      f_x = F_piece*dx;
		      f_y = F_piece*dy;
                        
		      forces->Fx[i] += f_x;
		      forces->Fy[i] += f_y;

		      
		      forces->Fx[j] += -f_x;
		      forces->Fy[j] += -f_y;
                                             
		      r_cm_a = fabs(-(p.BACTERIA_LENGTH - 1 - 2*a)*p.BALL_R);
		      r_cm_b = fabs(-(p.BACTERIA_LENGTH - 1 - 2*b)*p.BALL_R);
   
		      forces->Tau[i] += (f_y*cos(agents[i].th) 
					 - f_x*sin(agents[i].th) )*r_cm_a;

		      forces->Tau[j] += (-f_y*cos(agents[j].th) 
					 + f_x*sin(agents[j].th) )*r_cm_b;

                    }
                    
                }
            }

	  // printf("%d\t%lf\n", i, forces->Fx[i]);  
        }
    }
    
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

/***************************************************************************/
/*
  Clean out output files.
*/

void output_clean(struct Parameters p)
{
  FILE *fd = fopen("details.txt", "w");
  FILE *fx = fopen("cm_x_traj.txt", "w");
  FILE *fy = fopen("cm_y_traj.txt", "w");
  FILE *fth = fopen("th_traj.txt", "w");

  FILE *frob = fopen("xyt.txt", "w");

  // print a python dictionary of simulation details

  fprintf(fd, "{'r': %f, 'frame_lim': %f, 'run_time': %d, 'N': %d, 'L': %d}", 
	  p.BALL_R, p.SCREEN_W, p.RUN_TIME/p.SKIP, p.NUM_BACTERIA*p.BACTERIA_LENGTH, 
	  p.BACTERIA_LENGTH);
  
  fclose(frob);
  
  fclose(fth);
  fclose(fy);
  fclose(fx);
  fclose(fd);
}

/****************************************************************************/
/*
  Append data to text fields for animation with python
*/

void data_out(struct Parameters p, struct Agent *agents)
{   
  int i, j;

  FILE *fd = fopen("details.txt", "a");
  FILE *fx = fopen("cm_x_traj.txt", "a");
  FILE *fy = fopen("cm_y_traj.txt", "a");
  FILE *fth = fopen("th_traj.txt", "a");

  FILE *frob = fopen("xyt.txt", "a");
  
  for (i = 0; i < p.NUM_BACTERIA; i++)
    {
      for (j= 0; j < p.BACTERIA_LENGTH; j++)
	{
	  fprintf(fx, "%f , ", agents[i].balls[2*j]);
	  fprintf(fy, "%f , ", agents[i].balls[2*j+1]);

	  fprintf(frob, "%lf\t%lf\n", agents[i].balls[2*j], agents[i].balls[2*j+1]);
	}

      fprintf(fth, "%f , ", agents[i].th);
    }

  fprintf(frob,"\n\n");
 
  fclose(frob);
  fclose(fd);
  fclose(fx);
  fclose(fy);
  fclose(fth);
}  

void multiple_out(struct Parameters p, struct Agent *agents, int N)
{
  int i, j;
  
  char base[100], ind[100];

  FILE *fp;

  sprintf(ind,"%d", N);

  strcpy(base,"gnudat/time");

  strcat(base,ind);

  fp = fopen(base, "w");

  for (i = 0; i < p.NUM_BACTERIA; i++)
    {
      for (j= 0; j < p.BACTERIA_LENGTH; j++)
	{
	  fprintf(fp, "%lf\t%lf\n", agents[i].balls[2*j], agents[i].balls[2*j+1]);
	}
    }

  fclose(fp);
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
	  // data_out(p, agents);
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
