struct Parameters
{
  int RUN_TIME;
  int SKIP;
  int NUM_BACTERIA;
  int BACTERIA_LENGTH;
  
  int NUM_BOXES;
  int GRID_WIDTH; // num boxes in a row or column.  root of NUM_BOXES
  double BOX_WIDTH;
  
  int NPIL;
  double  PIL_SPAN;
  double PIL_LEN_MEAN;
  double PIL_LEN_STD;
  double PROB_EXTEND;
  double MOTOR_POWER;

  double K_STIFFNESS;
  double STATIC_FRICTION;
  double KINETIC_FRICTION;
  
  int UNIFORM;
  
  double SCREEN_W;
  double PLACEMENT_W;
  double BALL_R;
  double GAMMA;     //friction
  double E;         //energy scale

  int GRID;         // use grid (1) or standard MD (0).  Note Grid is nevessary for pilli to grab eachother...
  double DT;
};

struct Pillus
{
  //initial length, length, force, angle, motor power
  double L0, L, F, th, P;
  
  double s; // spring extension
  
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


/***********************************************************************************************
                                        organism.c
/***********************************************************************************************/
double * xy_position(struct Parameters p, int ID);
void make_agent(struct Parameters p, struct Agent * agents, int i, struct Box * grid, double x, double y, double th);
void make_colony(struct Parameters p, struct Agent *agents, long *idum, struct Box *grid);
void compute_rod(struct Parameters p, struct Agent * agents, int i);

void extend_pilli (struct Parameters p, long *idum, int i, struct Agent *agents);
void extend_pillus(struct Pillus * pil, struct Agent ag, long * idum, struct Parameters p);
void update_pilli(struct Agent * agents, int i, struct Parameters p);

/***********************************************************************************************
                                        forces.c
/***********************************************************************************************/
void friction(struct Agent * agents, int i, struct Parameters p);
void compute_pilli_forces(struct Agent * agents, int i, struct Parameters p);
void compute_forces(struct Parameters p, struct Agent *agents, double dt);
void verlet(double fx, double fy, double tau, struct Agent * agents, int i, double dt);

/***********************************************************************************************
                                        grid.c
/***********************************************************************************************/
void assign_grid_box(double box_width, int grid_width, int i, int j, struct Agent * agents, struct Box * grid);
void assign_grid_boxes(struct Parameters p, struct Agent * agents, int i, struct Box * grid);
void update_grid_position(struct Parameters p, int i, struct Agent *agents, struct Box * grid);

void compute_neighbours(struct Parameters p, int * neighbours, int grid_i);
void compute_forces_grid(struct Parameters p, struct Agent *agents, struct Box * grid, double dt);
/***********************************************************************************************
                                        main.c
/***********************************************************************************************/
void step(struct Parameters p, long *idum, int i, struct Agent *agents, double dt, int t);

void evolution(struct Parameters p, long *idum, 
	       struct Agent *agents, char * path, struct Box *grid);

/***********************************************************************************************
                                        misc.c
/***********************************************************************************************/
double compute_new_angle(double dx, double dy);
double min_sep(struct Parameters p, double a, double b);

double get_energy(struct Parameters p, struct Agent * agents, int norm);
void normalize_velocities(struct Parameters p, struct Agent * agents, double T, double V);



/***********************************************************************************************
                                        io.c
/***********************************************************************************************/
struct Parameters load_params(struct Parameters p);
void output_clean();
void mk_dirs(char * path);
void data_out(struct Parameters p, struct Agent *agents, char * path);
void multiple_out(struct Parameters, struct Agent *agents, int N, char * path);
void itoa(int n, char *s);

/***********************************************************************************************
                                        ran1.c
/***********************************************************************************************/
float ran1( long *idum);

