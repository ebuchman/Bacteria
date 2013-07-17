#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

/*****************************************************************************/
// methods

double wlc(double s, double L, double xi);
void compute_pilli_forces(struct Agent * agents, int i, struct Parameters p);
void update_pilli(struct Agent * agents, int i, struct Parameters p);
void extend_pilli (struct Parameters p, long *idum, int i, struct Agent *agents);
void extend_pillus(struct Pillus * pil, struct Agent ag, long *idum, struct Parameters p);

/*****************************************************************************/

void compute_pilli_forces(struct Agent * agents, int i, struct Parameters p)
{

  double fx, fy;

  int j;
  
  double s, L, th, T;

  fx = 0; fy = 0;

  for (j = 0; j < agents[i].Npil; j++)
  {
    s = agents[i].pillae[j].s;
    L = agents[i].pillae[j].L;
    th = agents[i].pillae[j].th;
  
    if (L > 0 && s > 0)
    {
	  T = wlc(s, L, p.XI);
      
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

void extend_pillus(struct Pillus * pil, struct Agent ag, long *idum, struct Parameters p)
{
  double r, t, dx, dy;
  double th_from_r, th;
  
  
  //uniform centred at mean with length 2*std
  pil->L0 = ran1(idum)*2*ag.pil_len_std + (ag.pil_len_mean - ag.pil_len_std);
  pil->L = pil->L0;
  pil->s = 0;
  
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

  pil->x = pbc(p, pil->x);
  pil->y = pbc(p, pil->y);

}

/*****************************************************************************/

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

/*****************************************************************************/

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
        pil->s += R - L; //extend sping
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
        pil->s += (pil->P / p.STATIC_FRICTION)*p.DT;
      }
      

      pil->th = compute_new_angle(dxp, dyp);
      
      
      int snapped;
      
      
      if (pil->s > 0.5*pil->L0) 
      { 
        //printf("snapped! x, L0: %f, %f\n\n", pil->s, pil->L0);
        pil->L = 0.0; // pillus snaps     
        snapped = 1;
      }
      
      if (pil->L <= 0.000000001)
      {
        if (snapped != 1) {printf("a pillus retracted completely\n"); }//exit(0);}
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
