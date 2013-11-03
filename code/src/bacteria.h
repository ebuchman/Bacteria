struct Parameters
{
  int RUN_TIME;
  int SKIP;
  int NUM_BACTERIA;
  int BACTERIA_LENGTH;
  
  int NUM_BOXES;
  int GRID_WIDTH; // num boxes in a row or column.
  double BOX_WIDTH;
  
  double SELF_FORCE;
  
  int NPIL;
  double  PIL_SPAN;
  double PIL_LEN_MEAN;
  double PIL_LEN_SD;
  double PROB_EXTEND;
  double MOTOR_POWER;

  double K_STIFFNESS;
  double STATIC_FRICTION;
  double KINETIC_FRICTION;
  
  int UNIFORM;
  
  double XI; 
  
  double SCREEN_W;
  double PLACEMENT_W;
  double BALL_R;
  double GAMMA;     //friction
  double E;         //energy scale

  int GRID;         // use grid (1) or standard MD (0).  Note Grid is nevessary for pilli to grab eachother
  int ATTACH;
  
  double DT;
};

struct Pillus
{
  //initial length, length, Tension, angle, motor power
  double L0, L, T, th, P;
  
  double s; // spring extension
  
  double x, y; //position of anchor

  int attached; // anchored to another agent?
  int agent_num; // which agent?
  int ball_num; // which ball?

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
  
  int * box_num; // location in the grid for each ball (single index)
  int agent_num;
  
  double iFx, iFy, iTau; // interaction forces
  double pFx, pFy, pTau; // pilli forces
  double fFx, fFy, fTau; // friction forces
    
    
  int t;
};

// each box in the grid has a pointer.  if there's a ball in that box, pointer points to agent and ball_num indexes ball.
struct Box
{
  int grid_pos; // position in grid .. necesssary?
  
  int occupied;
  
  int agent_num;
  int ball_num;
};


/* if i move this to misc.c, code will seg fault ... */
double pbc(struct Parameters p, double x);

/* if i move any of these, code will run weirdly ... */
double compute_new_angle(double dx, double dy);
double min_sep(struct Parameters p, double a, double b);
double pbc_th(double th);
void pbc_position(struct Agent * agents, int i, struct Parameters p);

float ran1( long *idum);

