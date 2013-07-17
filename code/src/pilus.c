#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

/*****************************************************************************/
// methods

double wlc(double s, double L, double xi);
//  void compute_pilli_forces(struct Agent * agents, int i, struct Parameters p);
void update_pilli(struct Agent * agents, int i, struct Parameters p);

/*****************************************************************************/

void compute_pilli_forces(struct Agent * agents, int i, struct Parameters p)
{

  double fx, fy;

  int j;
  
  double s, L, th, T;

  fx = 0; fy = 0;

  for (j = 0; j < p.NPIL; j++)
  {
    s = agents[i].pillae[j].s;
    L = agents[i].pillae[j].L;
    th = agents[i].pillae[j].th;
    //printf("i, j, T, L, s: %d, %d, %f, %f, %f\n", i, j, T, L, s);
  
    if (L > 0 && s > 0)
    {
	  T = wlc(s, L, p.XI);
      printf("i, j, T: %d, %d, %f\n", i, j, T);
      agents[i].pillae[j].T = T;
      
      fx += T*cos(th);
      fy += T*sin(th);
    }
          
    agents[i].pFx = fx;
    agents[i].pFy = fy;
    agents[i].pTau = p.BALL_R*p.BACTERIA_LENGTH*(-fx*sin(agents[i].th) + fy*cos(agents[i].th));
  }
}

/*****************************************************************************/

void extend_pilli (struct Parameters p, long *idum, struct Pillus *pil, double th0, double cmx, double cmy)
{
  int n;
  double r, t, dx, dy;
  double th;
  
  double eps = 1.0E-12;

  for (n = 0; n < p.NPIL; n++)
  {
    //printf("L: %f\n", pil[n].L);
    if (pil[n].L < eps)
    {

      if (ran1(idum) < p.PROB_EXTEND)
      {
        //uniform centred at mean with length 2*std
        pil[n].L0 = p.PIL_LEN_MEAN + (ran1(idum)*2*p.PIL_LEN_SD - p.PIL_LEN_SD);
        pil[n].L = pil[n].L0;
        pil[n].s = 0;
        
        //uniform around 0 with length pil_span
        th = th0 + 0.5*(2.0*ran1(idum)*p.PIL_SPAN - p.PIL_SPAN);
        th = pbc_th(th);
        
        pil[n].th = th;
              
        // anchor (x,y) : cm + to_end_of_rod + pilus extension  
        dx = pil[n].L*cos(th);
        dy = pil[n].L*sin(th);
        
        pil[n].x = cmx + p.BALL_R*p.BACTERIA_LENGTH*cos(th0) + dx;
        pil[n].y = cmy + p.BALL_R*p.BACTERIA_LENGTH*sin(th0) + dy;

        pil[n].x = pbc(p, pil[n].x);
        pil[n].y = pbc(p, pil[n].y);
        
        pil[n].P = p.MOTOR_POWER;
      }
    }
  }
}

/*****************************************************************************/

void update_pilli(struct Agent * agents, int i, struct Parameters p)
{
  double rod_end_x, rod_end_y;
  double dxp, dyp, R, L;
  double dU, dxs;
  double vx, vy;
  struct Pillus * pil;
  int j;
  
  double eps = 1.0E-12;

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
        pil->s += R - L; //extend sping
      }
      else //if (R < L)
      {
        pil->L += R - L; //retract pillus
        //printf("\n retract, R, L, dx, dy: %f, %f, %f, %f, %f\n", pil->L, R, L, dx, dy);
      }
      

      vx = fabs(agents[i].vx);
      vy = fabs(agents[i].vy);
      //printf("vx and vy: %f, %f\n", vx, vy);
      if (vx < eps && vy < eps)
      {
        printf("no motion\n");
        if (pil-> T < eps)
          pil->s += 0.1;
        else
          pil->s += (pil->P / pil->T)*p.DT;
      }

      pil->th = compute_new_angle(dxp, dyp);
      
      int snapped;
      
      
      if (pil->s > 0.3*pil->L) 
      { 
        //printf("snapped! x, L0: %f, %f\n\n", pil->s, pil->L0);
        pil->L = 0.0; // pillus snaps     
        snapped = 1;
      }
      
      if (pil->L <= eps)
      {
        pil->L = 0.0;
        pil->s = 0.0;
        pil->T = 0.0;   
      }

    }
  }
}

/*******************************************************************************/

double wlc(double s, double L, double xi)
{
  double f, x;

  double kT = 1.0;

  x = s/L;

  f = kT/xi*( 0.25/(1.0-x)/(1.0-x) - 0.25 + x);

  return f;
}
