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

void compute_forces_grid(struct Parameters p, struct Agent *agents, struct Box * grid, double dt)
{
  double L = p.BALL_R*2;
  
  
  int i, j, a; // iterators. i for agent, j for neighbour, a for ball
  
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
    agents[i].iFx = 0.0;
    agents[i].iFy = 0.0;
    agents[i].iTau = 0.0;
  }
  
  /* compute the new forces using the grid */
  
  /* for each agents, for each ball, consult grid position, check surrounding boxes for balls, compute forces */
  
  for (i = 0; i < p.NUM_BACTERIA; i++)
  {
      //printf("\n\n");

    for (a = 0; a < p.BACTERIA_LENGTH; a++)
    {
      grid_i = agents[i].box_num[a]; // grid index of the a'th ball
      compute_neighbours(p, grid_neighbours, grid_i);

      for (j = 0; j < 48; j ++)
      {
        this_box = grid_neighbours[j];
        this_agent = grid[this_box].agent_num;
        this_ball = grid[this_box].ball_num;
        
        // if grid box is occupied with a different agent
        if (grid[this_box].occupied == 1 && this_agent != i)
        {
          dx = min_sep(p, agents[i].balls[a*2], agents[this_agent].balls[this_ball*2]);
		  dy = min_sep(p, agents[i].balls[a*2 + 1], agents[this_agent].balls[this_ball*2 + 1]);
          
          //printf("i, a, b: %d, %d,%d \t dx, dy: %f, %f, ", i, a, this_ball, dx, dy);
          
		  r2 = dx*dx + dy*dy;
          
          if (r2 < pow(pow(2, 1./6)*L, 2)) // if close enough
          {
            F_piece = (48*p.E)*(pow(L, 12)*pow(r2, -6))/r2;
            
            if (F_piece > p.BALL_R/dt) // cap it for stability ...
              F_piece = p.BALL_R/dt; //pow(dt,2);
            
            f_x = F_piece*dx / 2.;
            f_y = F_piece*dy / 2.;
            
            //printf("fx, fy: %f, %f\n", f_x, f_y);
            
            agents[i].iFx += f_x;
            agents[i].iFy += f_y;
            
            agents[this_agent].iFx += -f_x;
            agents[this_agent].iFy += -f_y;
            
            r_cm_a = fabs(-(p.BACTERIA_LENGTH - 1 - 2*a)*p.BALL_R);
            r_cm_b = fabs(-(p.BACTERIA_LENGTH - 1 - 2*this_ball)*p.BALL_R);
            
            agents[i].iTau += ((f_y*cos(agents[i].th)
                               - f_x*sin(agents[i].th) )*r_cm_a)/2.;
            
            agents[this_agent].iTau += ((-f_y*cos(agents[this_agent].th)
                               + f_x*sin(agents[this_agent].th) )*r_cm_b)/2.;
          }
          //else printf("\n");
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
    //printf("%d, (%f, %f)\n", previous_box, agents[i].balls[2*j], agents[i].balls[2*j+1]);
    assign_grid_box(p.BOX_WIDTH, p.GRID_WIDTH, i, j, agents, grid);
    
    // if agent moved boxes
    if (agents[i].box_num[j] != previous_box)
    {
      // if agent moved but previous box still points to agent, unoccupy previous box
      if (grid[previous_box].agent_num == i)
        grid[previous_box].occupied = 0;
    }
      
  }//printf("\n\n");
}

/*****************************************************************************/
// assign boxes to agents and agents to boxes

void assign_grid_boxes(struct Parameters p, struct Agent * agents, int i, struct Box * grid)
{
  //printf("Box width: %f\n", p.BOX_WIDTH);

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

  x_box = (int) (x/box_width);
  y_box = (int) (y/box_width);
  
  box = y_box*grid_width + x_box;

  //printf("x, y, xbox, ybox, box: %f, %f, %d, %d, %d\n", x, y, x_box, y_box, box);

  // agent should point to grid and grid should point to agent (really they hold indices so we don't get tied up in pointers)
  agents[i].box_num[j] = box;

  grid[box].agent_num = i;
  grid[box].ball_num = j;
  
  grid[box].occupied = 1;
  
}

