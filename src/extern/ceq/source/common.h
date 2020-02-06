#ifndef common_h
#define common_h

extern const double TRACELIMIT;
extern const double tol;
extern const int attempts;

double update_limit_factor(double fac, double x, double dx);
void handle_trace_species_locking(double* a, double n, int nsp, int nel, double* ns, double* bi0, double* dlnns, int verbose);
int all_but_one_species_are_trace(int nsp, double* ns);
double constraint_errors(double* S, double* a,double* bi0,double* ns,int nsp,int nel,int neq,double* dlnns);
void composition_guess(double* a,double* M,double* X0,int nsp,int nel,double* ns,double* np,double* bi0);

#endif
