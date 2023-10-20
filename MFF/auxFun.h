#include "rel_euler.h"

void ciPP(struct PP eos, double vl,double rl,double vr,double rr, double *Pl, double *Pr, double *hl, double *hr);
void ciTASPP(struct TASPP eos, double vl,double rl,double vr,double rr, double *Pl, double *Pr, double *hl, double *hr);
void ciTAB(struct Tab eos, double vl,double rl,double vr,double rr, double *Pl, double *Pr, double *hl, double *hr);

void initialize(double*x,double**u,double **pri,double a,double b,int N, double rl, double vl, double Pl, double hl, double rr, double vr, double Pr, double hr);

void sumvecs(double *suma,double*a,double*b,int N);
void vectorbyn(double*sum,double s, double*a,int N);
double dot(double*a,double*b,int N);

void writePP(struct PP, double *x, double **u, double **pri, int N, char *name);
void writeTASPP(struct TASPP, double *x, double **u, double **pri, int N, char *name);
void writeTab(struct Tab, double *x, double **u, double **pri, int N, char* name);

void clean(void *ptr);
void cleanM(double **A, int N);

void displace(double **recons, int dim);