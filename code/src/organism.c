#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

/****************************************************************/

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

/*****************************************************************************/


void extend_pillus(struct Pillus * pil, int i, struct Agent ag, long *idum)
{
  //uniform centred at mean with length 2*std
  pil[i].L0 = ran1(idum)*2*ag.pil_len_std + (ag.pil_len_mean - ag.pil_len_std);

  //uniform around 0 with length pil_span
  pil[i].th = ran1(idum)*2*ag.pil_span - ag.pil_span;
  
}