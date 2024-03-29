#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <sys/stat.h>

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

  fscanf(fp, "%s %d", dummy, &(p.CONSERVATIVE) );
  fscanf(fp, "%s %lf", dummy, &(p.F_SELF) );  
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
  Make directory structure
*/

void mk_dirs(char * path)
{
  char pp[100];
  
  mkdir(path, S_IRWXU);
  
  sprintf(pp, "%s%s", path, "gnudat/");
  mkdir(pp, S_IRWXU);

  sprintf(pp, "%s%s", path, "pydat/");
  mkdir(pp, S_IRWXU);

}

/***************************************************************************/
/*
 Clean out output files.
 */

void output_clean(struct Parameters p, char * path)
{

  char pp[100];

  sprintf(pp, "%s%s", path, "pydat/details.txt");
  FILE *fd = fopen(pp, "w");
  
  sprintf(pp, "%s%s", path, "pydat/cm_x_traj.txt");
  FILE *fx = fopen(pp, "w");

  sprintf(pp, "%s%s", path, "pydat/cm_y_traj.txt");
  FILE *fy = fopen(pp, "w");
  
  sprintf(pp, "%s%s", path, "pydat/th_traj.txt");
  FILE *fth = fopen(pp, "w");
  
  sprintf(pp, "%s%s", path, "pydat/xyt.txt");
  FILE *frob = fopen(pp, "w");
  

  // print a python dictionary of simulation details
  
  fprintf(fd, "{'r': %f, 'frame_lim': %f, 'run_time': %d, 'skip': %d, 'N': %d, 'L': %d, 'dt': %f}",
          p.BALL_R, p.SCREEN_W, p.RUN_TIME, p.SKIP, p.NUM_BACTERIA,
          p.BACTERIA_LENGTH, p.DT);
  
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

void data_out(struct Parameters p, struct Agent *agents, char * path)
{
  int i, j;
  char pp[100];
  
  sprintf(pp, "%s%s", path, "pydat/details.txt");
  FILE *fd = fopen(pp, "a");
  
  sprintf(pp, "%s%s", path, "pydat/cm_x_traj.txt");
  FILE *fx = fopen(pp, "a");
  
  sprintf(pp, "%s%s", path, "pydat/cm_y_traj.txt");
  FILE *fy = fopen(pp, "a");
  
  sprintf(pp, "%s%s", path, "pydat/th_traj.txt");
  FILE *fth = fopen(pp, "a");
  
  sprintf(pp, "%s%s", path, "pydat/xyt.txt");
  FILE *frob = fopen(pp, "a");

  
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

void multiple_out(struct Parameters p, struct Agent *agents, int N, char * path)
{
  int i, j;
  
  char base[100], ind[100];
  char base2[100];
  
  FILE *fp;
  FILE *fp2;
  
  sprintf(ind,"%d", N);
  
  char pp[100];
  
  sprintf(pp, "%s%s", path, "gnudat/time");
  strcpy(base, pp);

  sprintf(pp, "%s%s", path, "gnudat/cm");
  strcpy(base2, pp);
  
  strcat(base,ind);
  strcat(base2,ind);
  
  fp = fopen(base, "w");
  fp2 = fopen(base2, "w");

  
  for (i = 0; i < p.NUM_BACTERIA; i++)
  {
  fprintf(fp2, "%f\t%f\n", agents[i].cm_x, agents[i].cm_y);
    
    
    for (j= 0; j < p.BACTERIA_LENGTH; j++)
	{
      fprintf(fp, "%f\t%f\n", agents[i].balls[2*j], agents[i].balls[2*j+1]);
      //printf("%f\t%f\n", agents[i].balls[2*j], agents[i].balls[2*j+1]);
     
	}
  }
  
  fclose(fp);
  fclose(fp2);
}
