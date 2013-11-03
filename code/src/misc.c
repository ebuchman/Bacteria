#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"




double get_energy(struct Parameters p, struct Agent * agents, int norm);
void normalize_velocities(struct Parameters p, struct Agent * agents, double T, double V);

/*************************************************/
// recenters atan to give [0, 2*PI]
double compute_new_angle(double dx, double dy)
{
  double phi; 
    
    
  if (dx == 0) 
  {
    if (dy > 0) phi = M_PI/2;
    else if (dy < 0) phi = M_PI*1.5;
    else phi = 0;
  }
  else
  {
    phi = atan(fabs(dy / dx));
  
    if (dx < 0 && dy > 0) phi = M_PI - phi;
    else if (dx < 0 && dy < 0) phi = M_PI + phi;
    else if (dx > 0 && dy < 0) phi = 2*M_PI - phi;
    else phi = phi;
  
  }
  return phi;

}

/*****************************************************************************/

//pbc min separation routine

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
//PBC

double pbc(struct Parameters p, double x)
{
  while (x<0.0)
    {
      x += p.SCREEN_W;
    }
     
  while (x >= p.SCREEN_W)
    {
      x -= p.SCREEN_W;
    }

  return x;
}

double pbc_th(double th)
{
  double pi2 = 2.0*M_PI;

  while (th<0.0)
    {
      th += pi2;
    }
     
  while (th >= pi2)
    {
      th -= pi2;
    }

  return th;
}

void pbc_position(struct Agent * agents, int i, struct Parameters p)
{
  agents[i].cm_x = pbc(p, agents[i].cm_x);
  agents[i].cm_y = pbc(p, agents[i].cm_y);
  agents[i].th = pbc_th(agents[i].th);
}

/*****************************************************************************/

double get_energy(struct Parameters p, struct Agent * agents, int norm)
{
  int i, j, a, b, N;
  double T, V, E; //kinetic, potential, total
  double dx, dy, r2;
  double L = p.BALL_R*2;
  
  N = p.NUM_BACTERIA;
  T = 0; V = 0; E = 0;
  
  
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
          
          printf("%f, ", r2);
          
          V += (4*p.E)*(pow(L, 12)*pow(r2, -6) -
                        pow(L, 6)*pow(r2, -3));
          
        }
      }
    }
    T += 0.5*(agents[i].vx*agents[i].vx + agents[i].vy*agents[i].vy);
  }
  
  if (norm ==1)
  {
    normalize_velocities(p, agents, T, V);
  }
  
  printf("\nV: %f, \tT: %f\n", V, T);
  
  E = T + V;
  
  return E;
}

void normalize_velocities(struct Parameters p, struct Agent * agents, double T, double V)
{
  int i, N;
  N = p.NUM_BACTERIA;
  double norm_factor;
  
  norm_factor = sqrt((p.E)/T);
  
  
  for (i=0 ; i < N; i++)
  {
    agents[i].vx = agents[i].vx * norm_factor;
    agents[i].vy = agents[i].vy * norm_factor;
  }
  
}