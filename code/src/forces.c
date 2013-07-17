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

void compute_pilli_forces(struct Agent * agents, int i, struct Parameters p)
{

  double this_fx, this_fy, this_tau;
  double x;
  int j;
  struct Pillus * pil;
  
  double s;

  agents[i].pFx = 0;
  agents[i].pFy = 0;
  agents[i].pTau = 0;

  for (j = 0; j < agents[i].Npil; j++)
  {
    pil = &agents[i].pillae[j];
    
    x = pil->s;
    if (pil->L > 0 && x > 0)
    {
      pil->F = p.K_STIFFNESS*x;
        
            
      this_fx = pil->F*cos(pil->th);
      this_fy = pil->F*sin(pil->th);
      this_tau = p.BALL_R*p.BACTERIA_LENGTH*(-this_fx*sin(agents[i].th) + this_fy*cos(agents[i].th));
    }
    else
    {
      this_fx = 0;
      this_fy = 0;
      this_tau = 0;
    }
      
    agents[i].pFx += this_fx;
    agents[i].pFy += this_fy;
    agents[i].pTau += this_tau;
  }
}

/*****************************************************************************/

// return vector of forces and vector of torques for entire colony

void compute_forces(struct Parameters p, struct Agent *agents, double dt)
{
  double L = p.BALL_R*2;

  int i, j, a, b; // iterators
  
  double dx, dy, r2;
  double F_piece, f_x, f_y;
  double r_cm_a, r_cm_b; //cm displacement's of balls
  
  
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
		  
		  dx = min_sep(p, agents[i].balls[a*2], agents[j].balls[b*2]);
		  
		  dy = min_sep(p, agents[i].balls[a*2 + 1], agents[j].balls[b*2 + 1]);
          
		  r2 = dx*dx + dy*dy;
          
		  if (r2 < pow(pow(2, 1./6)*L, 2)) //if close enough ...
          {
            //F_piece = (48*p.E)*(pow(L, 12)*pow(r2, -6) - 0.5*pow(L, 6)*pow(r2, -3))/r2;
            F_piece = (48*p.E)*(pow(L, 12)*pow(r2, -6))/r2;
            
		    
            if (F_piece > p.BALL_R/dt) // cap it for stability ...
              F_piece = p.BALL_R/dt;
            
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
      }
    }
  }
}
