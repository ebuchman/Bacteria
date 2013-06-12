struct Parameters
{
  int RUN_TIME;
  int SKIP;
  int NUM_BACTERIA;
  int BACTERIA_LENGTH;
  int UNIFORM;
  
  double SCREEN_W;
  double BALL_R;
  double GAMMA;     //friction
  double E;         //energy scale

  double DT;
};

struct Agent
{
  int N;
  double ball_r;
  
  double cm_x, cm_y;
  double th, omega;
  double *balls;
  
  double v0, F_self;
  double vx, vy;
  double last_Fx, last_Fy, last_tau;
  
  int t;
};

/* Convenient for returning from compute forces */

struct Forces
{
  double *Fx;
  double *Fy;
  double *Tau;
};

void evolution(struct Parameters p, long *idum, struct Forces *forces, 
	       struct Agent *agents);

void compute_forces(struct Parameters p, struct Forces *forces, 
		    struct Agent *agents, double dt);

void make_colony(struct Parameters p, struct Agent *agents, long *idum);

void make_colony_uniform(struct Parameters p, struct Agent *agents);

float ran1( long *idum);

void compute_rod(struct Parameters p, double *balls, double cm_x, double cm_y,
		 double r, double th, int N);


void step(struct Parameters p, long *idum, int i, struct Agent *agents, 
	  struct Forces *Fint, double dt);

double min_sep(struct Parameters p, double a, double b);

void output_clean();

void data_out(struct Parameters p, struct Agent *agents);

void multiple_out(struct Parameters, struct Agent *agents, int N);

void itoa(int n, char *s);
