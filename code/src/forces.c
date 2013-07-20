#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

/******************************************************************************/
void friction(struct Agent * agents, int i, struct Parameters p)
{
  double net_f, net_fx, net_fy;
  double fFx, fFy;
  double v, vx, vy;
  double epsilon = 1.0E-12;
  
  vx = agents[i].vx;
  vy = agents[i].vy;
  v = sqrt(vx*vx + vy*vy);
  
  net_fx = agents[i].iFx + agents[i].pFx;
  net_fy = agents[i].iFy + agents[i].pFy;
  net_f = sqrt(net_fx*net_fx + net_fy*net_fy);

  // static: if net force is greater than friction, otherwise, forces are 0
  //kinetic: if greater than friction, subtract friction.  otherwise, set to zero.  add velocity-dependent friction

  if (fabs(v) < epsilon)
  {      
      if (net_f < p.STATIC_FRICTION)
      {
        fFx = - net_fx;
        fFy = - net_fy;
      }
      else
      {
        fFx = - p.STATIC_FRICTION*net_fx/net_f;
        fFy = - p.STATIC_FRICTION*net_fy/net_f;
      }
  }
  else
  {
      fFx = - p.KINETIC_FRICTION*vx/v;
      fFy = - p.KINETIC_FRICTION*vy/v;
      
      /* velocity dependent friction */
      fFx -= p.GAMMA*vx;
      fFy -= p.GAMMA*vy;
  }
  
  agents[i].fFx = fFx;
  agents[i].fFy = fFy;
  /* We should calculate a frictional drag on the rotation */
}

/*****************************************************************************/
// compute interaction forces for all agents

void compute_forces(struct Parameters p, struct Agent *agents)
{
  int i, j, a, b; // iterators
  
  /* zero the forces */
  // would this be better at the end of step, cut down on an extra N loop?"
  for (i = 0; i < p.NUM_BACTERIA; i++)
  {
    agents[i].iFx = 0.0;
    agents[i].iFy = 0.0;
    agents[i].iTau = 0.0;
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
          force_core(agents, i, j, a, b, p);
        }
      }
    }
  }
}

/*****************************************************************************/
// compute forces for single interaction pair
// useful because we have two force routines, with and without the grid.  the force core is identical for both

void force_core(struct Agent *agents, int i, int j, int a, int b, struct Parameters p)
{
  double dx, dy, r2;
  double F_piece, f_x, f_y;
  double r_cm_a, r_cm_b; //cm displacement's of balls
  
  double L = p.BALL_R*2;
  
  dx = min_sep(p, agents[i].balls[a*2], agents[j].balls[b*2]);
  
  dy = min_sep(p, agents[i].balls[a*2 + 1], agents[j].balls[b*2 + 1]);
  
  r2 = dx*dx + dy*dy;
  
  if (r2 < pow(pow(2, 1./6)*L, 2)) //if close enough ...
  {
    //F_piece = (48*p.E)*(pow(L, 12)*pow(r2, -6) - 0.5*pow(L, 6)*pow(r2, -3))/r2; // LJ with attraction
    F_piece = (48*p.E)*(pow(L, 12)*pow(r2, -6))/r2; // LJ with just repulsion
    
    if (p.GRID == 1) F_piece = F_piece / 2; // since we will be double counting
    
    if (F_piece > p.BALL_R/p.DT) // cap it for stability ...
      F_piece = p.BALL_R/p.DT;
    
    f_x = F_piece*dx;
    f_y = F_piece*dy;
    
    agents[i].iFx += f_x;
    agents[i].iFy += f_y;
    
    agents[j].iFx += -f_x;
    agents[j].iFy += -f_y;
    
    r_cm_a = fabs(-(p.BACTERIA_LENGTH - 1 - 2*a)*p.BALL_R);
    r_cm_b = fabs(-(p.BACTERIA_LENGTH - 1 - 2*b)*p.BALL_R);
    
    agents[i].iTau += (f_y*cos(agents[i].th)
                       - f_x*sin(agents[i].th) )*r_cm_a;
    
    agents[j].iTau += (-f_y*cos(agents[j].th)
                       + f_x*sin(agents[j].th) )*r_cm_b;
  }
}
