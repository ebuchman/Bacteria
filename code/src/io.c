#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

/*****************************************************************************/

struct Parameters load_params(struct Parameters p)
{
  FILE *fp;
  
  char dummy[100];
  
  /* Read in parameters from input file */
  
  fp = fopen("input","r");
  
  fscanf(fp,"%s", dummy);
  
  fscanf(fp, "%s %d", dummy, &(p.RUN_TIME) );
  fscanf(fp, "%s %d", dummy, &(p.SKIP) );
  fscanf(fp, "%s %d", dummy, &(p.NUM_BACTERIA) );
  fscanf(fp, "%s %d", dummy, &(p.BACTERIA_LENGTH) );
  
  fscanf(fp, "%s %d", dummy, &(p.NPIL) );
  fscanf(fp, "%s %lf", dummy, &(p.PIL_SPAN) );
  fscanf(fp, "%s %lf", dummy, &(p.PIL_LEN_MEAN) );
  fscanf(fp, "%s %lf", dummy, &(p.PIL_LEN_STD) );
  
  fscanf(fp, "%s %lf", dummy, &(p.K_STIFFNESS) );
  fscanf(fp, "%s %lf", dummy, &(p.F_FRICTION) );
  
  fscanf(fp, "%s %d", dummy, &(p.UNIFORM) );
  
  fscanf(fp, "%s %lf", dummy, &(p.SCREEN_W) );
  fscanf(fp, "%s %lf", dummy, &(p.BALL_R) );
  fscanf(fp, "%s %lf", dummy, &(p.GAMMA) );
  fscanf(fp, "%s %lf", dummy, &(p.E) );
  
  fscanf(fp, "%s %lf", dummy, &(p.DT) );
  
  fclose(fp);
    
  return p;
  
}

/***************************************************************************/
/*
 Clean out output files.
 */

void output_clean(struct Parameters p)
{
  FILE *fd = fopen("pydat/details.txt", "w");
  FILE *fx = fopen("pydat/cm_x_traj.txt", "w");
  FILE *fy = fopen("pydat/cm_y_traj.txt", "w");
  FILE *fth = fopen("pydat/th_traj.txt", "w");
  
  FILE *frob = fopen("pydat/xyt.txt", "w");
  
  // print a python dictionary of simulation details
  
  fprintf(fd, "{'r': %f, 'frame_lim': %f, 'run_time': %d, 'N': %d, 'L': %d}",
          p.BALL_R, p.SCREEN_W, p.RUN_TIME/p.SKIP, p.NUM_BACTERIA*p.BACTERIA_LENGTH,
          p.BACTERIA_LENGTH);
  
  fclose(frob);
  
  fclose(fth);
  fclose(fy);
  fclose(fx);
  fclose(fd);
}

/****************************************************************************/
/*
 Append data to text fields for animation with python
 */

void data_out(struct Parameters p, struct Agent *agents)
{
  int i, j;
  
  FILE *fd = fopen("pydat/details.txt", "a");
  FILE *fx = fopen("pydat/cm_x_traj.txt", "a");
  FILE *fy = fopen("pydat/cm_y_traj.txt", "a");
  FILE *fth = fopen("pydat/th_traj.txt", "a");
  
  FILE *frob = fopen("pydat/xyt.txt", "a");
  
  for (i = 0; i < p.NUM_BACTERIA; i++)
  {
    for (j= 0; j < p.BACTERIA_LENGTH; j++)
	{
	  fprintf(fx, "%f , ", agents[i].balls[2*j]);
	  fprintf(fy, "%f , ", agents[i].balls[2*j+1]);
      
	  fprintf(frob, "%lf\t%lf\n", agents[i].balls[2*j], agents[i].balls[2*j+1]);
	}
    
    fprintf(fth, "%f , ", agents[i].th);
  }
  
  fprintf(frob,"\n\n");
  
  fclose(frob);
  fclose(fd);
  fclose(fx);
  fclose(fy);
  fclose(fth);
}

void multiple_out(struct Parameters p, struct Agent *agents, int N)
{
  int i, j;
  
  char base[100], ind[100];
  
  FILE *fp;
  
  sprintf(ind,"%d", N);
  
  
  strcpy(base,"gnudat/time");
  
  strcat(base,ind);
  
  fp = fopen(base, "w");
  
  for (i = 0; i < p.NUM_BACTERIA; i++)
  {
    for (j= 0; j < p.BACTERIA_LENGTH; j++)
	{
      fprintf(fp, "%f\t%f\n", agents[i].balls[2*j], agents[i].balls[2*j+1]);
      //printf("%f\t%f\n", agents[i].balls[2*j], agents[i].balls[2*j+1]);
      
	}
  }
  
  fclose(fp);
}
