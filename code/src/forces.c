#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

/*****************************************************************************/

void compute_pilli_forces(struct pilForces * forces, struct Agent * agents, int i, struct Parameters p)
{

  double this_fx, this_fy, this_tau;
  double x;
  int j;
  struct Pillus * pil;
  
  double s;

  for (j = 0; j < agents[i].Npil; j++)
  {
    pil = &agents[i].pillae[j];
    
    x = pil->x_ext;
    //printf ("%f, ", x);
    if (pil->L > 0 && x > 0)
    {
      pil->F = p.K_STIFFNESS*x;
        
            
      this_fx = pil->F*cos(pil->th);
      this_fy = pil->F*sin(pil->th);
      this_tau = p.BALL_R*p.BACTERIA_LENGTH*(-this_fx*sin(agents[i].th) + this_fy*cos(agents[i].th));
      
      
      forces->Fx += this_fx;
      forces->Fy += this_fy;
      forces->Tau += this_tau;
      
      
      //forces->Fx[i] += this_fx;
      //forces->Fy[i] += this_fy;
      //forces->Tau[i] += this_tau;
      
      
      
      
      //printf("fx, fy, tau: %f, %f, %f\n", this_fx, this_fy, this_tau);

    }
  }//printf("\n");

    
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
    }
  }
  
}
