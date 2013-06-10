struct Agent{
    int N;
    double ball_r;

    double cm_x, cm_y;
    double th, omega;
    double v_r, v_t;
    double * balls;
    
    double F_self;
    double last_F_r, last_F_t, last_tau;
    
    int t;
};

/* Convenient for returning from compute forces */
struct Forces{
    double * F_r;
    double * F_t;
    double * Tau;

};


struct Agent * make_agent(int N, double x, double y, double th); 
struct Agent ** make_colony(); 
double * compute_rod(double cm_x, double cm_y, double r, double th, int N); 
void step(struct Agent * agent, double F_r, double F_t, double tau, double dt);
struct Forces * compute_forces(struct Agent ** agents, double dt); 