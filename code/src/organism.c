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
  double V, W, A, N, L, D, offset, wx, wy;
  int nx, ny;
  double x, y, dx, dy;
  
  W = p.PLACEMENT_W;
  V = W*W;
  N = p.NUM_BACTERIA;
  A = V/N;
  
  D = p.BALL_R*2;
  L = D*p.BACTERIA_LENGTH;
  
  offset = (-(L+D) + sqrt((L-D)*(L-D) + 4*A))/4.; //solution to quadratic optimization/placement problem (first attempt)
  
  wx = D + 2*offset;
  wy = L + 2*offset;
  
  nx = (int) W / wx;
  ny = (int) W / wy;
  
  // now that we have an approximate fit, we force everything to fit and no overlaps.
  while (nx*ny < N)
  {
    nx+=2;
    ny+=1;
  }
  
  //recompute offset and everything
  // N is effectively bigger, making A smaller - but there will probably be left over room not filled
  N = nx*ny;
  A = V/N;

  offset = (-(L+D) + sqrt((L-D)*(L-D) + 4*A))/4.; //solution to quadratic optimization/placement problem
  wx = D + 2*offset;
  wy = L + 2*offset;
  
  x = (int) (ID%nx)*wx + (offset + p.BALL_R);
  y = (int) (ID/nx)*wy + (offset + L/2 );
  
  //centre
  dx = (p.SCREEN_W - W)/2.;
  dy = (p.SCREEN_W - W)/2.;
  
  x+= dx;
  y+=dy;
  
  xy[0] = x;
  xy[1] = y;
  
  return xy;
}


/****************************************************************/
// Make single agent
// This is useful if we want bacteria to be able to divide

void make_agent(struct Parameters p, struct Agent * agents, int i, struct Box * grid, double x, double y, double th)
{
  int j;
    
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
      
  agents[i].ball_r = p.BALL_R; // this is dumb...
  
  agents[i].iFx = 0;
  agents[i].iFy = 0;
  agents[i].iTau = 0;

  agents[i].pFx = 0;
  agents[i].pFy = 0;
  agents[i].pTau = 0;
  
  agents[i].fFx = 0;
  agents[i].fFy = 0;
  agents[i].fTau = 0;
  
  agents[i].last_Fx = 0;
  agents[i].last_Fy = 0;
  
  agents[i].last_tau = 0;
  
  agents[i].t = 0;
  
  compute_rod(p, agents, i);
  
  if (p.GRID == 1)
    assign_grid_boxes(p, agents, i, grid);
      
  for (j=0; j < agents[i].Npil; j++)
  {
    agents[i].pillae[j].P = p.MOTOR_POWER;
  
  }  

}

/****************************************************************/
// Make colony
void make_colony(struct Parameters p, struct Agent *agents, long *idum, struct Box *grid)
{
  int i;
  double x, y, th;
  
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
    make_agent(p, agents, i, grid, x, y, th);
  }
}

/*****************************************************************************/
/*
 extend an agent from his cm in direction of theta
 */

//void compute_rod(struct Parameters p, double *balls, double cm_x, double cm_y, double r, double th, int N)
void compute_rod(struct Parameters p, struct Agent * agents, int n)
{
  int i, N;
  double x, y, cm_x, cm_y, r, th;
  
  cm_x = agents[n].cm_x;
  cm_y = agents[n].cm_y;
  
  th = agents[n].th;
  r = p.BALL_R;
  N = agents[n].N;
  
  for (i=0; i<N; i++)
  {
    
    x = cm_x - (N - 1 - 2*i)*r*cos(th);
    y = cm_y - (N - 1 - 2*i)*r*sin(th);
    
    x = pbc(p, x);
    y = pbc(p, y);
    
    agents[n].balls[i*2] = x;
    agents[n].balls[i*2+1] = y;
  }
}

/*****************************************************************************/




// Reverse on bounce back (concept):

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
