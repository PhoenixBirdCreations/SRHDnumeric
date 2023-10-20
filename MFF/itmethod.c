#include "itmethod.h"

void recoverPP(struct PP eos, double **u, double **pri, int N){
 double x,D,S,tau;
 int j,iter; 
 double w,rho;
 for(j=1;j<N+1;j++){
    D=u[j][0]; S=u[j][1]; tau=u[j][2];
    x=0; iter=0;    
    while(x>=0 && iter<20){ 
      if(S>(tau+x+D)) {printf("theorem violated: D %f S %f Tau %f, index %d\n",D,S,tau,j); exit(1);}
      w=1/sqrt(1-S*S/((tau+x+D)*(tau+x+D)));
      rho=D/w;
      x=pressurePP(rho,eos);
      iter++;
    }
    if(x<=0){
     printf("D=%f S=%f tau=%f , x=%f\n",D,S,tau,x);
     printf("tuple (D,S,tau) was not physic, in position %d. Recovery failed.\n",j);
     exit(1);
    }
    else{
     pri[j][1]=S/(tau+D+x); pri[j][0]=D*sqrt(1-pri[j][1]*pri[j][1]); pri[j][2]=x; 
    }
 }  
}

void recoverTASPP(struct TASPP eos, double **u, double **pri, int N){
 double x,D,S,tau;
 int j,iter; 
 double w,rho;
 for(j=1;j<N+1;j++){
    D=u[j][0]; S=u[j][1]; tau=u[j][2];
    x=0; iter=0;    
    while(x>=0 && iter<20){ 
      if(S>(tau+x+D)) {printf("theorem violated: D %f S %f Tau %f, index %d\n",D,S,tau,j); exit(1);}
      w=1/sqrt(1-S*S/((tau+x+D)*(tau+x+D)));
      rho=D/w;
      x=pressureTASPP(rho,eos);
      iter++;
    }
    if(x<=0){
     printf("D=%f S=%f tau=%f , x=%f\n",D,S,tau,x);
     printf("tuple (D,S,tau) was not physic, in position %d. Recovery failed.\n",j);
     exit(1);
    }
    else{
     pri[j][1]=S/(tau+D+x); pri[j][0]=D*sqrt(1-pri[j][1]*pri[j][1]); pri[j][2]=x; 
    }
 }  
}

void recoverTab(struct Tab eos, double **u, double **pri, int N){
 double x,D,S,tau;
 int j,iter; 
 double w,rho;
 for(j=1;j<N+1;j++){
    D=u[j][0]; S=u[j][1]; tau=u[j][2];
    x=0; iter=0;    
    while(x>=0 && iter<20){ 
      if(S>(tau+x+D)) {printf("theorem violated: D %f S %f Tau %f, index %d\n",D,S,tau,j); exit(1);}
      w=1/sqrt(1-S*S/((tau+x+D)*(tau+x+D)));
      rho=D/w;
      x=pressureTab(rho,eos);
      iter++;
    }
    if(x<=0){
     printf("D=%f S=%f tau=%f , x=%f\n",D,S,tau,x);
     printf("tuple (D,S,tau) was not physic, in position %d. Recovery failed.\n",j);
     exit(1);
    }
    else{
     pri[j][1]=S/(tau+D+x); pri[j][0]=D*sqrt(1-pri[j][1]*pri[j][1]); pri[j][2]=x; 
    }
 }  
}