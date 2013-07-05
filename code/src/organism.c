#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

/****************************************************************/
// method to set uniform initial density
double * xy_position(struct Parameters p, int ID)
{
  double * xy = malloc(sizeof(double)*2);
  double V, A, N, L, D, offset, wx, wy;
  int nx, ny;
  double x, y;
  
  V = p.SCREEN_W*p.SCREEN_W;
  N = p.NUM_BACTERIA;
  A = V/N;
  
  D = p.BALL_R*2;
  L = D*p.BACTERIA_LENGTH;
  
  offset = (-(L+D) + sqrt((L-D)*(L-D) + 4*A))/4.;
  
  wx = D + 2*offset;
  wy = L + 2*offset;
  
  nx = (int) p.SCREEN_W / wx; //round up
  ny = (int) p.SCREEN_W / wy;
  
  if (nx*ny < N)
  {
    nx+=1;
    ny+=1;
  }
  
  x = (ID%nx)*wx + (offset + p.BALL_R);
  y = (int) (ID/nx)*wy + (offset + L/2 );
  
  xy[0] = x;
  xy[1] = y;
  
  return xy;
}

/****************************************************************/

void make_colony(struct Parameters p, struct Agent *agents, long *idum)
{
  double  x, y, th;
  int i, j;
  
  for (i=0; i < p.NUM_BACTERIA; i++)
  {
    
    if (p.UNIFORM == 0)
    {
      x = p.SCREEN_W*ran1(idum);
      y = p.SCREEN_W*ran1(idum);
      th = 2.0*M_PI*ran1(idum);
    }
    else
    {
      double * xy;
      xy = xy_position(p, i);
      x = xy[0];
      y = xy[1];
      th = M_PI/2;
    }
    
    agents[i].N = p.BACTERIA_LENGTH;
    agents[i].Npil = p.NPIL;
    agents[i].pil_span = p.PIL_SPAN;
    agents[i].pil_len_mean = p.PIL_LEN_MEAN;
    agents[i].pil_len_std = p.PIL_LEN_STD;
    
    agents[i].cm_x = x;
    agents[i].cm_y = y;
    
    agents[i].th = th;
    
    agents[i].omega = 0;
    agents[i].vx = 0;
    agents[i].vy = 0;
        
    agents[i].ball_r = p.BALL_R;
    
    agents[i].last_Fx = 0;
    agents[i].last_Fy = 0;
    
    agents[i].last_tau = 0;
    
    agents[i].t = 0;
    
    compute_rod(p, agents[i].balls, agents[i].cm_x, agents[i].cm_y,
                agents[i].ball_r, agents[i].th, agents[i].N);
    
    
    for (j=0; j < agents[i].Npil; j++)
    {
      agents[i].pillae[j].P = p.MOTOR_POWER;
    
    }
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


void extend_pillus(struct Pillus * pil, int i, struct Agent ag, long *idum, struct Parameters p)
{
  double r, t, dx, dy;
  
  //uniform centred at mean with length 2*std
  pil[i].L0 = ran1(idum)*2*ag.pil_len_std + (ag.pil_len_mean - ag.pil_len_std);
  pil[i].L = pil[i].L0;
  pil[i].x_ext = 0;
  
  //uniform around 0 with length pil_span
  pil[i].th = ran1(idum)*ag.pil_span - ag.pil_span/2;
  
  // anchor (x,y) : cm + to_end_of_rod + pilus extension  
  r = pil[i].L * cos(pil[i].th);
  t = pil[i].L * sin(pil[i].th);
  
  dx = r*cos(ag.th) - t*sin(ag.th);
  dy = r*sin(ag.th) + t*cos(ag.th);
  
  pil[i].x = ag.cm_x + p.BALL_R*p.BACTERIA_LENGTH*cos(ag.th) + dx;
  pil[i].y = ag.cm_y + p.BALL_R*p.BACTERIA_LENGTH*sin(ag.th) + dy;
  
  // pbc
  pil[i].x = fmod(pil[i].x, p.SCREEN_W);
  pil[i].y = fmod(pil[i].y, p.SCREEN_W);
  
  printf("LO: %f, th: %f, x:, %f, y: %f\n", pil[i].L, pil[i].th, pil[i].x, pil[i].y);
  
}