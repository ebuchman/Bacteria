#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

/*****************************************************************************/

void compute_pilli_forces(struct Forces * forces, struct Agent * agents, int i, struct Parameters p)
{

  double f_x=0; double f_y=0; double tau=0;
  double f_r=0; double f_t=0;
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
      if (pil->F > p.F_FRICTION) pil->F = p.F_FRICTION;
        
      f_r += pil->F*cos(pil->th);
      f_t += pil->F*sin(pil->th);
      tau += p.BALL_R*p.BACTERIA_LENGTH*pil->F*sin(pil->th);
    }
  }//printf("\n");
  forces->Fx[i] += f_r*cos(agents[i].th) - f_t*sin(agents[i].th);
  forces->Fy[i] += f_r*sin(agents[i].th) + f_t*cos(agents[i].th) ;
  forces->Tau[i] += tau;
  
  printf("Pillus: x: %f, th: %f, F: %f, f_r: %f, f_t: %f\n", x, pil->th, pil->F, f_r, f_t);
      
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
