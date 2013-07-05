struct Parameters
{
  int RUN_TIME;
  int SKIP;
  int NUM_BACTERIA;
  int BACTERIA_LENGTH;
    
  int NPIL;
  double  PIL_SPAN;
  double PIL_LEN_MEAN;
  double PIL_LEN_STD;
  double PROB_EXTEND;
  double MOTOR_POWER;

  double K_STIFFNESS;
  double F_FRICTION;

  int UNIFORM;
  
  double SCREEN_W;
  double BALL_R;
  double GAMMA;     //friction
  double E;         //energy scale

  double DT;
};

struct Pillus
{
  //initial length, length, force, angle, motor power
  double L0, L, F, th, P;
  
  double x_ext; // spring extension
  
  double x, y; //position of anchor

};

struct Agent
{
  int N;
  double ball_r;
  
  double cm_x, cm_y;
  double th, omega;
  double *balls;
  
  double vx, vy;
  double last_Fx, last_Fy, last_tau;
    
  int Npil;
  double pil_span, pil_len_mean, pil_len_std;
  struct Pillus *pillae;
  
    
  int t;
};

/* Convenient for returning from compute forces */

struct Forces
{
  double *Fx;
  double *Fy;
  double *Tau;
};


double * xy_position(struct Parameters p, int ID);

double compute_new_angle(double dx, double dy, double th1);

double get_energy(struct Parameters p, struct Agent * agents, int norm);

void normalize_velocities(struct Parameters p, struct Agent * agents, double T, double V);

struct Parameters load_params(struct Parameters p);

void evolution(struct Parameters p, long *idum, struct Forces *forces, 
	       struct Agent *agents, char * path);

void extend_pillus(struct Pillus * pi, int i, struct Agent ag, long * idum, struct Parameters p);

void compute_pilli_forces(struct Forces * forces, struct Agent * agents, int i, struct Parameters p);

void compute_forces(struct Parameters p, struct Forces *forces, 
		    struct Agent *agents, double dt);

void make_colony(struct Parameters p, struct Agent *agents, long *idum);

void make_colony_uniform(struct Parameters p, struct Agent *agents);

float ran1( long *idum);

void compute_rod(struct Parameters p, double *balls, double cm_x, double cm_y,
		 double r, double th, int N);


void step(struct Parameters p, long *idum, int i, struct Agent *agents, 
	  struct Forces *Fint, double dt, int t);

double min_sep(struct Parameters p, double a, double b);

void output_clean();

void mk_dirs(char * path);

void data_out(struct Parameters p, struct Agent *agents, char * path);

void multiple_out(struct Parameters, struct Agent *agents, int N, char * path);

void itoa(int n, char *s);
