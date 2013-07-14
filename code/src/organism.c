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

void make_colony(struct Parameters p, struct Agent *agents, long *idum, struct Box *grid)
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
        
    agents[i].ball_r = p.BALL_R; // this is dumb...
    
    agents[i].last_Fx = 0;
    agents[i].last_Fy = 0;
    
    agents[i].last_tau = 0;
    
    agents[i].t = 0;
    
    //compute_rod(p, agents[i].balls, agents[i].cm_x, agents[i].cm_y, agents[i].ball_r, agents[i].th, agents[i].N);
    compute_rod(p, agents, i);
    
    assign_grid_boxes(p, agents, i, grid);
        
    for (j=0; j < agents[i].Npil; j++)
    {
      agents[i].pillae[j].P = p.MOTOR_POWER;
    
    }
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
    
    // pbc
    
    x = fabs(fmod(x + 3.0*p.SCREEN_W, p.SCREEN_W));
    y = fabs(fmod(y + 3.0*p.SCREEN_W, p.SCREEN_W));
    
    agents[n].balls[i*2] = x;
    agents[n].balls[i*2+1] = y;
  }
}

/*****************************************************************************/


void extend_pillus(struct Pillus * pil, struct Agent ag, long *idum, struct Parameters p)
{
  double r, t, dx, dy;
  double th_from_r, th;
  
  
  //uniform centred at mean with length 2*std
  pil->L0 = ran1(idum)*2*ag.pil_len_std + (ag.pil_len_mean - ag.pil_len_std);
  pil->L = pil->L0;
  pil->x_ext = 0;
  
  //uniform around 0 with length pil_span
  th_from_r = ran1(idum)*ag.pil_span - ag.pil_span/2; //angle from radius of agent
  th = ag.th + th_from_r;  //angle in xy

  if (th < 0) th = 2*M_PI + th;
  else if (th > 2*M_PI) th = th - 2*M_PI;
  
  pil->th = th;
    
        
  // anchor (x,y) : cm + to_end_of_rod + pilus extension  
  dx = pil->L*cos(th);
  dy = pil->L*sin(th);
  
  pil->x = ag.cm_x + p.BALL_R*p.BACTERIA_LENGTH*cos(ag.th) + dx;
  pil->y = ag.cm_y + p.BALL_R*p.BACTERIA_LENGTH*sin(ag.th) + dy;

    
  // pbc
  pil->x = fmod(pil->x, p.SCREEN_W);
  pil->y = fmod(pil->y, p.SCREEN_W);
  
  //printf("LO: %f, th: %f, x:, %f, y: %f\n", pil[i].L, pil[i].th, pil[i].x, pil[i].y);
  
}

void extend_pilli (struct Parameters p, long *idum, int i, struct Agent *agents)
{
  int j;

  for (j = 0; j < agents[i].Npil; j++)
  {
    if (agents[i].pillae[j].L <= 0.00000001)
    {

      if (ran1(idum) < p.PROB_EXTEND)
      {
        
        extend_pillus(&agents[i].pillae[j], agents[i], idum, p);

      }
    }
  }
}

void update_pilli(struct Agent * agents, int i, struct Parameters p)
{
  double rod_end_x, rod_end_y;
  double dxp, dyp, R, L;
  double dU, dxs;
  double vx, vy;
  struct Pillus * pil;
  int j;

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
      
      if (R >= L)
      {
        pil->x_ext += R - L; //extend sping
      }
      else //if (R < L)
      {
        pil->L += R - L; //retract pillus
        //printf("\n retract, R, L, dx, dy: %f, %f, %f, %f, %f\n", pil->L, R, L, dx, dy);
      }
      

      vx = fabs(agents[i].vx);
      vy = fabs(agents[i].vy);

      if (vx < 0.000001 && vy < 0.000001)
      {
        //printf("no motion\n");
        pil->x_ext += (pil->P / p.STATIC_FRICTION)*p.DT;
      }
      

      pil->th = compute_new_angle(dxp, dyp);
      
      
      int snapped;
      
      
      if (pil->x_ext > 0.5*pil->L0) 
      { 
        //printf("snapped! x, L0: %f, %f\n\n", pil->x_ext, pil->L0);
        pil->L = 0.0; // pillus snaps     
        snapped = 1;
      }
      
      if (pil->L <= 0.000000001)
      {
        if (snapped != 1) {printf("a pillus retracted completely\n"); }//exit(0);}
        pil->L = 0.0;
        pil->x_ext = 0.0;
        pil->F = 0.0;   
      }

    }
  }
}


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
