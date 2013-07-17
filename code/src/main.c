#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

int main()
{
  int i;

  long *idum;

  struct Agent *agents;
  
  struct Box *grid; // an array of boxes, re-presenting grid positions

  struct Parameters p;

  char path[100];

  idum = (long *)malloc(sizeof(long));
  *idum = (long) 271; // consistent
  //*idum = -(long) time(NULL); // random 

  /* Load Parameters */

  p = load_params(p);
  
  /* Build Directory Structure */
  sprintf(path, "data/N%d_T%d/", p.NUM_BACTERIA, p.RUN_TIME);
  printf("here: %s\n", path);
  mk_dirs(path); // make directories for data
    
  /*Initialize grid*/
  if (p.GRID == 1)
  {
    grid = (struct Box *)malloc(sizeof(struct Box)*p.NUM_BOXES);
    for (i=0; i<p.NUM_BOXES; i++)
    {
      grid[i].occupied = 0;
      grid[i].grid_pos = i;
    }
  }
  
  /* Initialize agents */
  
  agents = (struct Agent *)malloc(sizeof(struct Agent)*p.NUM_BACTERIA);

  // malloc for grid pointers, balls xy's, and pillae
  for(i=0; i<p.NUM_BACTERIA; i++)
    {
      agents[i].box_num = (int *)malloc(sizeof(int)*p.BACTERIA_LENGTH);
      agents[i].balls = (double *)malloc(p.BACTERIA_LENGTH*2*sizeof(double));
      agents[i].pillae = (struct Pillus *)malloc(p.NPIL*sizeof(struct Pillus));
    }    

  /* Initialize colony variables */
  make_colony(p, agents, idum, grid); 

  /*** RUN SIMULATION ***/
  printf("evolution\n");
  evolution(p, idum, agents, path, grid);

  printf("# toodles!\n");

  free(idum);
    
  free(agents);

  return 0;

}