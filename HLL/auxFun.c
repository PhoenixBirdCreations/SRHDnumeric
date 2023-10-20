#include "auxFun.h"

void ciTASPP(struct TASPP eos, double vl,double rl,double vr,double rr, double *Pl, double *Pr, double *hl, double *hr){
  *Pl=pressureTASPP(rl,eos); *Pr=pressureTASPP(rr,eos);
  *hl=hTASPP(rl,eos); *hr=hTASPP(rr,eos);
}
void ciPP(struct PP eos, double vl,double rl,double vr,double rr, double *Pl, double *Pr, double *hl, double *hr){
  *Pl=pressurePP(rl,eos); *Pr=pressurePP(rr,eos);
  *hl=hPP(rl,eos); *hr=hPP(rr,eos);
}
void ciTAB(struct Tab eos, double vl,double rl,double vr,double rr, double *Pl, double *Pr, double *hl, double *hr){
  *Pl=pressureTab(rl,eos); *Pr=pressureTab(rr,eos);
  *hl=hTab(rl,eos); *hr=hTab(rr,eos);
}

void initialize(double*x,double**u,double **pri,double a,double b,int N, double rl, double vl, double Pl, double hl, double rr, double vr, double Pr, double hr){
  double wl, wr; wl=W(vl); wr=W(vr);
  double c1l=rl*wl, c1r=rr*wr;
  double c2l=rl*hl*wl*wl*vl,    c3l=rl*hl*wl*wl-Pl-c1l;
  double c2r=rr*hr*wr*wr*vr,    c3r=rr*hr*wr*wr-Pr-c1r;
  double dx=(b-a)/(N-1);
  int i;
  for(i=0;i<N/2;i++){
    x[i]=a+i*dx;
    x[N-1-i]=b-i*dx;

    u[i][0]=c1l; u[i][1]=c2l; u[i][2]=c3l;
    pri[i][0]=rl; pri[i][1]=vl; pri[i][2]=Pl;

    u[N+5-i][0]=c1r; u[N+5-i][1]=c2r; u[N+5-i][2]=c3r;
    pri[N+5-i][0]=rr; pri[N+5-i][1]=vr; pri[N+5-i][2]=Pr;
  }
  
  x[N/2]=a+(N/2)*dx;
  u[N/2][0]=c1l;   u[N/2][1]=c2l;   u[N/2][2]=c3l;
  u[N/2+1][0]=c1l; u[N/2+1][1]=c2l; u[N/2+1][2]=c3l;
  u[N/2+2][0]=c1l; u[N/2+2][1]=c2l; u[N/2+2][2]=c3l;
  u[N/2+3][0]=c1r; u[N/2+3][1]=c2r; u[N/2+3][2]=c3r;
  u[N/2+4][0]=c1r; u[N/2+4][1]=c2r; u[N/2+4][2]=c3r;
  u[N/2+5][0]=c1r; u[N/2+5][1]=c2r; u[N/2+5][2]=c3r;
  u[N/2+6][0]=c1r; u[N/2+6][1]=c2r; u[N/2+6][2]=c3r;

  pri[N/2][0]=rl;   pri[N/2][1]=vl;   pri[N/2][2]=Pl; 
  pri[N/2+1][0]=rl; pri[N/2+1][1]=vl; pri[N/2+1][2]=Pl; 
  pri[N/2+2][0]=rl; pri[N/2+2][1]=vl; pri[N/2+2][2]=Pl; 
  pri[N/2+3][0]=rr; pri[N/2+3][1]=vr; pri[N/2+3][2]=Pr; 
  pri[N/2+4][0]=rr; pri[N/2+4][1]=vr; pri[N/2+4][2]=Pr; 
  pri[N/2+5][0]=rr; pri[N/2+5][1]=vr; pri[N/2+5][2]=Pr; 
  pri[N/2+6][0]=rr; pri[N/2+6][1]=vr; pri[N/2+6][2]=Pr; 
}

void sumvecs(double *suma,double*a,double*b,int N){
  int i;
  for (i=0;i<N;i++){
    suma[i]=a[i]+b[i];
  }
}  
void vectorbyn(double*sum,double s, double*a,int N){
 int i;
 for(i=0;i<N;i++){
  sum[i]=a[i]*s; 
 }
} 
double dot(double*a,double*b,int N){
  int i;double sum=0;
  for (i=0;i<N;i++){
    sum+=a[i]*b[i];
  }
  return sum;
}

void writePP(struct PP eos,double*x,double**u,double **pri,int N,char*name){
  FILE * fp;
  fp = fopen (name, "w+");
  int i;
  fprintf(fp, "#1-position  2-density   3-velocity    4-pressure    5-fund_derivative   6-nonlinearityfactor\n");
  for (i=0;i<N;i++){
    fprintf(fp, "%e %.10e %.10e %.10e %.10e %.10e\n", x[i], pri[i+1][0], pri[i+1][1], pri[i+1][2], GPP(pri[i+1][0],eos), nonlfactorPP(pri[i+1][0],eos));
  }
  fclose(fp);
}
void writeTASPP(struct TASPP eos,double*x,double**u,double **pri,int N,char*name){
  FILE * fp;
  fp = fopen (name, "w+");
  int i;
  fprintf(fp, "#1-position  2-density   3-velocity    4-pressure    5-fund_derivative   6-nonlinearityfactor\n"); 
  for (i=0;i<N;i++){
    fprintf(fp, "%e %.10e %.10e %.10e %.10e %.10e\n", x[i], pri[i+1][0], pri[i+1][1], pri[i+1][2], GTASPP(pri[i+1][0],eos), nonlfactorTASPP(pri[i+1][0],eos));
  }
  fclose(fp);
}
void writeTab(struct Tab eos,double*x,double**u,double **pri,int N,char*name){
  FILE * fp;
  fp = fopen (name, "w+");
  int i;
  fprintf(fp, "#1-position  2-density   3-velocity    4-pressure    5-fund_derivative   6-nonlinearityfactor\n"); 
  for (i=0;i<N;i++){
    fprintf(fp,"%e %.10e %.10e %.10e %.10e %.10e\n", x[i], pri[i+1][0], pri[i+1][1], pri[i+1][2], GTab(pri[i+1][0],eos), nonlfactorTab(pri[i+1][0],eos));
  }
  fclose(fp);
}

void clean(void *ptr){
 free(ptr);ptr=NULL; 
}
void cleanM(double **A, int N){
  int i;
  for(i=0;i<N;i++)
    free(A[i]);
  free(A);
}

void displace(double **recons, int dim){
  int i;
  for (i=0; i<4; i++)
    memcpy(recons[i], recons[i+2], dim*sizeof(double));
}