#ifndef ceq_h
#define ceq_h

int pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1,
       int verbose);

int rhou(double rho,double u,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1,
         double* T,int verbose);

int ps(double pt,double st,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
         double* X1, double* T, int verbose);

double get_u(double T, double* X, int nsp, double* lewis, double* M);

double get_h(double T, double* X, int nsp, double* lewis, double* M);

double get_cp(double T, double* X, int nsp, double* lewis, double* M);

double get_s0(double T, double* X, int nsp, double* lewis, double* M);

double get_s(double T, double p, double* X, int nsp, double* lewis, double* M);

int batch_pt(int N, double* p,double* T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
             double* X1, int verbose);

int batch_rhou(int N, double* rho,double* u,double* X0,int nsp,int nel,double* lewis,double* M,
               double* a, double* X1, double* T, int verbose);

int batch_u(int N, double* T, double* X, int nsp, double* lewis, double* M, double* u);

#endif
