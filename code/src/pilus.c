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

void compute_pilli_forces(struct Agent * agents, struct Parameters p)
{

  double fx, fy, tau;

  int i, j;
  
  double s, L, th, T;

  // zero the forces
  for (i=0; i < p.NUM_BACTERIA; i++)
  {
    agents[i].pFx = 0;
    agents[i].pFy = 0;
    agents[i].pTau = 0;
  }
  
  for (i = 0; i < p.NUM_BACTERIA; i++)
  {
    fx = 0; fy = 0; tau = 0;
    for (j = 0; j < p.NPIL; j++)
    {
      s = agents[i].pillae[j].s;
      L = agents[i].pillae[j].L;
      th = agents[i].pillae[j].th;
      //printf("i, j, T, L, s: %d, %d, %f, %f, %f\n", i, j, T, L, s);
    
      if (L > 0 && s > 0)
      {
        //T = p.K_STIFFNESS*s;
        T = wlc(s, L, p.XI);
        //printf("i, j, T: %d, %d, %f\n", i, j, T);
        agents[i].pillae[j].T = T;
        
        fx = T*cos(th);
        fy = T*sin(th);
        tau = p.BALL_R*p.BACTERIA_LENGTH*(-fx*sin(agents[i].th) + fy*cos(agents[i].th));
      }
      
      agents[i].pFx += fx;
      agents[i].pFy += fy;
      agents[i].pTau += tau;   
      
      // add forces from pilli that attach to other bacteria
      if (p.ATTACH == 1 && agents[i].pillae[j].attached == 1)
      {
        int agent_num, ball_num;
        
        agent_num = agents[i].pillae[j].agent_num;
        ball_num = agents[i].pillae[j].ball_num;
        
        agents[agent_num].pFx += -fx;
        agents[agent_num].pFy += -fy;
        agents[agent_num].pTau += p.BALL_R*(ball_num - (p.BACTERIA_LENGTH - 1)/2)*(fx*sin(agents[agent_num].th) - fy*cos(agents[agent_num].th)) ;

      }
    }
  }
}

/*****************************************************************************/

void extend_pilli (struct Parameters p, long *idum, struct Box *grid, struct Agent *agents, int i, struct Pillus *pil, double th0, double cmx, double cmy)
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
        pil[n].P = p.MOTOR_POWER;

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
        
        // check if pillus is anchored to another bacteria (requires grid)
        // reset pillus anchor to centre of ball in grid box
        if (p.GRID == 1 && p.ATTACH == 1)
        {
          int box, ball_num, agent_num;
          double rod_end_x, rod_end_y;
          double dxp, dyp;
          
          box = box_from_xy(pil[n].x, pil[n].y, p.BOX_WIDTH, p.GRID_WIDTH);
          
          if (grid[box].occupied == 1)
          {
            agent_num = grid[box].agent_num;
            ball_num = grid[box].ball_num;
            
            //printf("attachmen!\n");
            
            // pillus includes agent its attached to
            pil[n].attached = 1;
            pil[n].agent_num = agent_num;
            pil[n].ball_num = ball_num;
            
            //reset xy, length, etc. 
            pil[n].x = agents[agent_num].balls[2*ball_num];
            pil[n].y = agents[agent_num].balls[2*ball_num + 1];
            
            rod_end_x = agents[i].cm_x + p.BALL_R*p.BACTERIA_LENGTH*cos(agents[i].th);
            rod_end_y = agents[i].cm_y + p.BALL_R*p.BACTERIA_LENGTH*sin(agents[i].th);
            
            rod_end_x = pbc(p, rod_end_x);
            rod_end_y = pbc(p, rod_end_y);

            dxp = min_sep(p, pil[n].x, rod_end_x);
            dyp = min_sep(p, pil[n].y, rod_end_y);
            
            pil[n].L0 = sqrt(dxp*dxp + dyp*dyp);
            pil[n].L = pil[n].L0;
            pil[n].th = compute_new_angle(dxp, dyp);
            
          }
        }
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
  double retract_v;
  struct Pillus * pil;
  int j;
  
  double eps = 1.0E-12;

  rod_end_x = agents[i].cm_x + p.BALL_R*p.BACTERIA_LENGTH*cos(agents[i].th);
  rod_end_y = agents[i].cm_y + p.BALL_R*p.BACTERIA_LENGTH*sin(agents[i].th);
  
  rod_end_x = pbc(p, rod_end_x);
  rod_end_y = pbc(p, rod_end_y);

  
  for (j=0; j < agents[i].Npil; j++)
  {
    
    pil = &agents[i].pillae[j];
    
    if (pil->L > eps)
    {
      // if attached to other agent, update anchor point to new cm of other agent's ball
      if (pil->attached == 1)
      {
        pil->x = agents[pil->agent_num].balls[2*pil->ball_num];
        pil->y = agents[pil->agent_num].balls[2*pil->ball_num + 1];
      }
      
      dxp = min_sep(p, pil->x, rod_end_x);
      dyp = min_sep(p, pil->y, rod_end_y);
      R = sqrt(dxp*dxp + dyp*dyp);
      L = pil->L;
      
      /* If the new position is further from the anchor point, we extend the spring;'
            If it is closer, we retract the pillus.  */
      
      if (R - L > eps)
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
        //printf("no motion\n");
        
        if (pil-> T < p.STATIC_FRICTION)
          retract_v = pil->P/p.STATIC_FRICTION;
        else
          retract_v = pil->P/pil->T;
        
        pil->s += retract_v*p.DT;
        
        
        
        
        //pil->s += (pil->P / p.STATIC_FRICTION)*p.DT;
      }

      pil->th = compute_new_angle(dxp, dyp);
      
      int snapped;
      
      
      if (pil->s > 0.8*pil->L) 
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
        pil->attached = 0;
      }
    }
    //  printf("T: %f\t", pil->T);
    
  }
  //printf("\n");
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
