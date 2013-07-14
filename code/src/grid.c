#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"


// 48 neighbours (7x7 grid around centre)
void compute_neighbours(struct Parameters p, int * neighbours, int grid_i)
{
  int i;
  int r, c, w; //row num, column num, grid width
   
  w = p.GRID_WIDTH;
  c = grid_i - r*w;
  
  //pbc'd
  
  //top row
  r = (int) (grid_i / w + 3)%w;
  for(i = 0; i<7; i++)
    neighbours[i] = r*w + (c + i - 3)%w;  

  //second row above centre
  r = (int) (grid_i / w + 2)%w;
  for(i = 0; i<7; i++)
    neighbours[i + 7] = r*w + (c + i - 3)%w;  
    
  //first row above centre
  r = (int) (grid_i / w + 1)%w;
  for(i = 0; i<7; i++)
    neighbours[i + 14] = r*w + (c + i - 3)%w;    
  
  // centre
  r = (int) grid_i / w;
  for (i=0; i<3; i++)
    neighbours[i + 21] = r*w + (c + i - 3)%w;
  for (i=0; i<3; i++)
    neighbours[i + 24] = r*w + (c + i + 1)%w;
      
  //first row below centre
  r = (int) (grid_i / w - 1)%w;
  for(i = 0; i<7; i++)
    neighbours[i + 27] = r*w + (c + i - 3)%w ;  

  //second row below centre
  r = (int) (grid_i / w - 2)%w;
  for(i = 0; i<7; i++)
    neighbours[i + 34] = r*w + (c + i - 3)%w ; 

  //bottom row
  r = (int) (grid_i / w - 3)%w;
  for(i = 0; i<7; i++)
    neighbours[i + 41] = r*w + (c + i - 3)%w; 

}


/*****************************************************************************/

// return vector of forces and vector of torques for entire colony
// computed using the grid (N time complexity instead of N^2)

void compute_forces_grid(struct Parameters p, struct Forces *forces,
                    struct Agent *agents, struct Box * grid, double dt)
{
  double L = p.BALL_R*2;
  
  
  int i, j, a; // iterators
  
  int grid_i; //grid index of agent from main loop
  int this_agent, this_ball, this_box; // index agent (and ball) and box for interactions (neighbours)
  
  int * grid_neighbours; // array of indices for neighbour boxes
  
  grid_neighbours = (int *)malloc(sizeof(int)*48);
  
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
  
  /* compute the new forces using the grid */
  
  /* for each agents, for each ball, consult grid position, check surrounding boxes for balls, compute forces */
  
  for (i = 0; i < p.NUM_BACTERIA; i++)
  {
    for (a = 0; a < p.BACTERIA_LENGTH; a++)
    {
      grid_i = agents[i].box_num[a]; // grid index of the a'th ball
      compute_neighbours(p, grid_neighbours, grid_i);

      for (j = 0; j < 48; j ++)
      {
        this_box = grid_neighbours[j];
        this_agent = grid[this_box].agent_num;
        this_ball = grid[this_box].ball_num;
        
        if (grid[this_box].occupied == 1 && this_agent != i)
        {
          dx = min_sep(p, agents[i].balls[a*2], agents[this_agent].balls[this_ball*2]);
		  dy = min_sep(p, agents[i].balls[a*2 + 1], agents[this_agent].balls[this_ball*2 + 1]);
          
		  r2 = dx*dx + dy*dy;
          
          F_piece = (48*p.E)*(pow(L, 12)*pow(r2, -6) -
                              0.5*pow(L, 6)*pow(r2, -3))/r2;
          
          if (F_piece > p.BALL_R/dt) // cap it for stability ...
            F_piece = p.BALL_R/dt; //pow(dt,2);
          
          f_x = F_piece*dx;
          f_y = F_piece*dy;
          
          forces->Fx[i] += f_x;
          forces->Fy[i] += f_y;
          
          forces->Fx[this_agent] += -f_x;
          forces->Fy[this_agent] += -f_y;
          
          r_cm_a = fabs(-(p.BACTERIA_LENGTH - 1 - 2*a)*p.BALL_R);
          r_cm_b = fabs(-(p.BACTERIA_LENGTH - 1 - 2*this_ball)*p.BALL_R);
          
          forces->Tau[i] += (f_y*cos(agents[i].th)
                             - f_x*sin(agents[i].th) )*r_cm_a;
          
          forces->Tau[this_agent] += (-f_y*cos(agents[this_agent].th)
                             + f_x*sin(agents[this_agent].th) )*r_cm_b;
        }
      }
    
    }
  }
  
  free(grid_neighbours);
}

/*****************************************************************************/
// update grid positions after a step

void update_grid_position(struct Parameters p, int i, struct Agent *agents, struct Box * grid)
{
  // use centre of mass of balls to compute grid position
  int j; 
  int previous_box; 
  
  
  for (j = 0; j < p.BACTERIA_LENGTH; j++)
  {
    previous_box = agents[i].box_num[j];
    assign_grid_box(p.BOX_WIDTH, p.GRID_WIDTH, i, j, agents, grid);
    
    // if agent moved boxes
    if (agents[i].box_num[j] != previous_box)
    {
      // if agent moved but previous box still points to agent, unoccupy previous box
      if (grid[previous_box].agent_num == i)
        grid[previous_box].occupied = 0;
    }
      
  }
}

/*****************************************************************************/
// assign boxes to agents and agents to boxes

void assign_grid_boxes(struct Parameters p, struct Agent * agents, int i, struct Box * grid)
{
  // use centre of mass of balls to compute grid position
  int j; 
    
  for (j = 0; j < p.BACTERIA_LENGTH; j++)
  {
    assign_grid_box(p.BOX_WIDTH, p.GRID_WIDTH, i, j, agents, grid);
  }
}

/*****************************************************************************/
// compute and assign grid position for single ball of single agent

void assign_grid_box(double box_width, int grid_width, int i, int j, struct Agent * agents, struct Box * grid)
{
  
  double x, y;
  int x_box, y_box, box;
  
  x = agents[i].balls[2*j];   // i is the agent num
  y = agents[i].balls[2*j + 1];   // j is the ball num

  x_box = (int) x/box_width;
  y_box = (int) y/box_width;
  
  box = y_box*grid_width + x_box;

  // agent should point to grid and grid should point to agent (really they hold indices so we don't get tied up in pointers)
  agents[i].box_num[j] = box;

  grid[box].agent_num = i;
  grid[box].ball_num = j;
  
  grid[box].occupied = 1;
}
