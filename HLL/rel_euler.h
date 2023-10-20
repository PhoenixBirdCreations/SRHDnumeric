#include "itmethod.h"

double W(double x);
void f(double *re,double*u, double *pri);

double eigenvaluePP(int k, double *pri, struct PP eos);
double eigenvalueTASPP(int k, double *pri, struct TASPP eos);
double eigenvalueTab(int k, double *pri, struct Tab eos);

void eigenvectorsPP(struct PP eos, double **r, double **l, double*pl, double*pr);
void eigenvectorsTASPP(struct TASPP eos, double **r, double **l, double*pl, double*pr);
void eigenvectorsTab(struct Tab eos, double **r, double **l, double*pl, double*pr);

