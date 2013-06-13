#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "bacteria.h"

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
