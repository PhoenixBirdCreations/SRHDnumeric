#include "auxFun.h"

int scalarsPP(struct PP eos, int allowUpwind, int method, double dx, double *phiplus, double* phimin, int j, double **u,double **pri, double**l,int arg, int N);
int scalarsTASPP(struct TASPP eos, int allowUpwind, int method, double dx, double *phiplus, double* phimin, int j, double **u,double **pri, double**l,int arg, int N);
int scalarsTab(struct Tab eos, int allowUpwind, int method, double dx, double *phiplus, double* phimin, int j, double **u,double **pri, double**l,int arg, int N);

void FHLL_PP(struct PP eos, double* result,int j, double **u, double** pri, double *reconL, double *reconR);
void FHLL_TASPP(struct TASPP eos, double* result,int j, double **u, double** pri, double *reconL, double *reconR);
void FHLL_TAB(struct Tab eos, double* result,int j, double **u, double** pri, double *reconL, double *reconR);

void integratePP(struct PP eos, int allowUpwind, int method, int N, double**in, double **out, double **pri, double dx, double cfl, double**l, double **r);
void integrateTASPP(struct TASPP eos, int allowUpwind, int method, int N, double**in, double **out, double **pri, double dx, double cfl, double**l, double **r);
void integrateTab(struct Tab eos, int allowUpwind, int method, int N, double**in, double **out, double **pri, double dx, double cfl, double**l, double **r);


void rk3(int usePP, int useTASPP, void *eosGeneric, double *x, int allowUpwind, int method, double tf, int N, double **in, double **out, double **pri, double dx, double dt, double cfl, double **l, double **r);