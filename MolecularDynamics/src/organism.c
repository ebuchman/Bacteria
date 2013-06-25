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
  
  if (p.CONSERVATIVE==1)
  {
    // hamiltonian.  No friction.  No driving.
    p.F_SELF = 0;
    p.F_FRICTION = 0;
    p.GAMMA = 0;
  }
  
  for (i=0; i < p.NUM_BACTERIA; i++)
  {
    agents[i].N = p.BACTERIA_LENGTH;

    if (p.UNIFORM == 0)
    {
      x = p.SCREEN_W*ran1(idum);
      y = p.SCREEN_W*ran1(idum);
    }
    else
    {
      double * xy;
      xy = xy_position(p, i);
      x = xy[0];
      y = xy[1];
      printf("%f, %f\n", x, y);
    }
    agents[i].cm_x = x;
    agents[i].cm_y = y;
    
    agents[i].omega = 0;
    
    if (p.CONSERVATIVE == 0)
    {
      agents[i].vx = 0;
      agents[i].vy = 0;
    }
    else
    {
      agents[i].vx = ran1(idum)*2 - 1; // velocities between +/- 1.  To be normalized
      agents[i].vy = ran1(idum)*2 - 1;
    }
    
    
    if (p.BACTERIA_LENGTH > 1)
    {
      if (p.UNIFORM == 0)
      {
        th = 2.0*M_PI*ran1(idum);
        agents[i].th = th;
      }
      else
        agents[i].th = M_PI/2; // everyone points up
    }
    else
      agents[i].th = atan(agents[i].vy/agents[i].vx);

    
    agents[i].ball_r = p.BALL_R;
    
    agents[i].last_Fx = 0;
    agents[i].last_Fy = 0;
    
    agents[i].last_tau = 0;
    
    agents[i].t = 0;
    
    compute_rod(p, agents[i].balls, agents[i].cm_x, agents[i].cm_y,
                agents[i].ball_r, agents[i].th, agents[i].N);
    
  }
  
  if (p.CONSERVATIVE != 0)
  {
    get_energy(p, agents, 1); // 1 tells it to normalize
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
