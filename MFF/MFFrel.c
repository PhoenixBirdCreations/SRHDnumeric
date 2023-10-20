#include "MFFrel.h"

#define SIGN(x) (((x)>0) ? 1 : -1)

int scalarsPP(struct PP eos, int allowUpwind, int method, double dx, double *phiplus, double* phimin, int j, double **u,double **pri, double**l,int arg, int N){ 
  int i,k;
  int dim=3;
  double *valuef, *empty, **m;
  if((valuef = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory valuef\n"); exit(41);}
  if((empty = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory 'empty'\n"); exit(41);}
  if ((m=(double**) malloc((2*method+1)*sizeof(double*)))==NULL) {printf("Error allocating memory for reconstruction\n"); exit(42);}
  for (i=0; i<(2*method+1); i++) if((m[i]=(double*) malloc(dim*sizeof(double)))==NULL) {printf("Error allocating memory\n"); exit(42);}
  

  if (allowUpwind && !arg){
    int signizq,signder;
    double l1i=eigenvaluePP(1,pri[j],eos), l2i=eigenvaluePP(2,pri[j],eos), l3i=eigenvaluePP(3,pri[j],eos), l1d=eigenvaluePP(1,pri[j+1],eos), l2d=eigenvaluePP(2,pri[j+1],eos), l3d=eigenvaluePP(3,pri[j+1],eos);
    signizq=SIGN(l1i)+SIGN(l2i)+SIGN(l3i);
    signder=SIGN(l1d)+SIGN(l2d)+SIGN(l3d);
    if(l1i*l1d>0 && l2i*l2d>0 && l3i*l3d>0 && (signizq==3 || signizq==-3) && (signder==3 || signder==-3)){ //just if possible
      double *phiplustemp, *phimintemp;
      if((phiplustemp = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory phipt\n"); exit(41);}
      if((phimintemp = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory phimt\n"); exit(41);}

      for(i=0;i<2*method+1;i++){
        for(k=0;k<dim;k++){
          f(valuef, u[j-method+i],pri[j-method+i]);
          m[i][k]=dot(l[k],valuef,dim);
        }
      }
      reconstruct(method,m,method,dim,dx,phiplustemp,empty);

      for(i=0;i<2*method+1;i++){
        for(k=0;k<dim;k++){
          f(valuef, u[j-method+i+1],pri[j-method+i+1]);
          m[i][k]=dot(l[dim+k],valuef,dim);
        }
      }
      reconstruct(method,m,method,dim,dx,empty,phimintemp);

      for (k=0;k<3;k++){
        if(eigenvaluePP(k+1,pri[j],eos)>0) {
          phiplus[k]=phiplustemp[k];
          phimin[k]=0.0;
        }
        else{
          phiplus[k]=0.0;
          phimin[k]=phimintemp[k];
        }
      }
      cleanM(m,(2*method+1)); free(valuef); free(phimintemp); free(phiplustemp); free(empty);
      return 0;
    }
  }
  
  double *alpha; double max1, max2, max3; int state, statem1;
  if((alpha = (double*)malloc((4+2*(method-1))*sizeof(double))) == NULL) {printf("Error allocating memory alpha\n"); exit(41);}
  for(i=0;i<4+2*(method-1);i++){
     state=-method+i; if (j+state>=N+5) statem1=state; else statem1=state+1;
     max1=fmax(fabs(eigenvaluePP(1,pri[j+state],eos)),fabs(eigenvaluePP(1,pri[j+statem1],eos)));
     max2=fmax(fabs(eigenvaluePP(2,pri[j+state],eos)),fabs(eigenvaluePP(2,pri[j+statem1],eos)));
     max3=fmax(fabs(eigenvaluePP(3,pri[j+state],eos)),fabs(eigenvaluePP(3,pri[j+statem1],eos)));
     max1=fmax(max1,max2);
     alpha[i]=fmax(max1,max3);
     if(arg) alpha[i]=sqrt(0.5*(1+alpha[i]*alpha[i]));
  }
  for(i=1;i<4+2*(method-1);i++) alpha[0]=fmax(alpha[0],alpha[i]);

  for(i=0;i<2*method+1;i++){
    for(k=0;k<dim;k++){
      f(valuef, u[j-method+i],pri[j-method+i]);
      m[i][k]=0.5*(dot(l[k],valuef,dim)+alpha[0]*dot(l[k],u[j-method+i],dim));
    }
  }
  reconstruct(method, m, method, dim, dx, phiplus, empty);

  for(i=0;i<2*method+1;i++){
    for(k=0;k<dim;k++){
      f(valuef, u[j-method+i+1],pri[j-method+i+1]);
      m[i][k]=0.5*(dot(l[dim+k],valuef,dim)-alpha[0]*dot(l[dim+k],u[j-method+i+1],dim));
    }
  }
  reconstruct(method, m, method, dim, dx, empty, phimin);

  cleanM(m, 2*method+1); free(alpha); free(empty); free(valuef);
  return 0;
  
}
int scalarsTASPP(struct TASPP eos, int allowUpwind, int method, double dx, double *phiplus, double* phimin, int j, double **u,double **pri, double**l,int arg, int N){ 
  int i,k;
  int dim=3;
  double *valuef, *empty, **m;
  if((valuef = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory valuef\n"); exit(41);}
  if((empty = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory 'empty'\n"); exit(41);}
  if ((m=(double**) malloc((2*method+1)*sizeof(double*)))==NULL) {printf("Error allocating memory for reconstruction\n"); exit(42);}
  for (i=0; i<(2*method+1); i++) if((m[i]=(double*) malloc(dim*sizeof(double)))==NULL) {printf("Error allocating memory\n"); exit(42);}
  

  if (allowUpwind && !arg){
    int signizq,signder;
    double l1i=eigenvalueTASPP(1,pri[j],eos), l2i=eigenvalueTASPP(2,pri[j],eos), l3i=eigenvalueTASPP(3,pri[j],eos), l1d=eigenvalueTASPP(1,pri[j+1],eos), l2d=eigenvalueTASPP(2,pri[j+1],eos), l3d=eigenvalueTASPP(3,pri[j+1],eos);
    signizq=SIGN(l1i)+SIGN(l2i)+SIGN(l3i);
    signder=SIGN(l1d)+SIGN(l2d)+SIGN(l3d);
    if(l1i*l1d>0 && l2i*l2d>0 && l3i*l3d>0 && (signizq==3 || signizq==-3) && (signder==3 || signder==-3)){ //just if possible
      double *phiplustemp, *phimintemp;
      if((phiplustemp = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory phipt\n"); exit(41);}
      if((phimintemp = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory phimt\n"); exit(41);}

      for(i=0;i<2*method+1;i++){
        for(k=0;k<dim;k++){
          f(valuef, u[j-method+i],pri[j-method+i]);
          m[i][k]=dot(l[k],valuef,dim);
        }
      }
      reconstruct(method,m,method,dim,dx,phiplustemp,empty);

      for(i=0;i<2*method+1;i++){
        for(k=0;k<dim;k++){
          f(valuef, u[j-method+i+1],pri[j-method+i+1]);
          m[i][k]=dot(l[dim+k],valuef,dim);
        }
      }
      reconstruct(method,m,method,dim,dx,empty,phimintemp);

      for (k=0;k<3;k++){
        if(eigenvalueTASPP(k+1,pri[j],eos)>0) {
          phiplus[k]=phiplustemp[k];
          phimin[k]=0.0;
        }
        else{
          phiplus[k]=0.0;
          phimin[k]=phimintemp[k];
        }
      }
      cleanM(m,(2*method+1)); free(valuef); free(phimintemp); free(phiplustemp); free(empty);
      return 0;
    }
  }
  
  double *alpha; double max1, max2, max3; int state, statem1;
  if((alpha = (double*)malloc((4+2*(method-1))*sizeof(double))) == NULL) {printf("Error allocating memory alpha\n"); exit(41);}
  for(i=0;i<4+2*(method-1);i++){
     state=-method+i; if (j+state>=N+5) statem1=state; else statem1=state+1;
     max1=fmax(fabs(eigenvalueTASPP(1,pri[j+state],eos)),fabs(eigenvalueTASPP(1,pri[j+statem1],eos)));
     max2=fmax(fabs(eigenvalueTASPP(2,pri[j+state],eos)),fabs(eigenvalueTASPP(2,pri[j+statem1],eos)));
     max3=fmax(fabs(eigenvalueTASPP(3,pri[j+state],eos)),fabs(eigenvalueTASPP(3,pri[j+statem1],eos)));
     max1=fmax(max1,max2);
     alpha[i]=fmax(max1,max3);
     if(arg) alpha[i]=sqrt(0.5*(1+alpha[i]*alpha[i]));
  }
  for(i=1;i<4+2*(method-1);i++) alpha[0]=fmax(alpha[0],alpha[i]);

  for(i=0;i<2*method+1;i++){
    for(k=0;k<dim;k++){
      f(valuef, u[j-method+i],pri[j-method+i]);
      m[i][k]=0.5*(dot(l[k],valuef,dim)+alpha[0]*dot(l[k],u[j-method+i],dim));
    }
  }
  reconstruct(method, m, method, dim, dx, phiplus, empty);

  for(i=0;i<2*method+1;i++){
    for(k=0;k<dim;k++){
      f(valuef, u[j-method+i+1],pri[j-method+i+1]);
      m[i][k]=0.5*(dot(l[dim+k],valuef,dim)-alpha[0]*dot(l[dim+k],u[j-method+i+1],dim));
    }
  }
  reconstruct(method, m, method, dim, dx, empty, phimin);

  cleanM(m, 2*method+1); free(alpha); free(empty); free(valuef);
  return 0;
  
}
int scalarsTab(struct Tab eos, int allowUpwind, int method, double dx, double *phiplus, double* phimin, int j, double **u,double **pri, double**l,int arg, int N){ 
  int i,k;
  int dim=3;
  double *valuef, *empty, **m;
  if((valuef = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory valuef\n"); exit(41);}
  if((empty = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory 'empty'\n"); exit(41);}
  if ((m=(double**) malloc((2*method+1)*sizeof(double*)))==NULL) {printf("Error allocating memory for reconstruction\n"); exit(42);}
  for (i=0; i<(2*method+1); i++) if((m[i]=(double*) malloc(dim*sizeof(double)))==NULL) {printf("Error allocating memory\n"); exit(42);}
  

  if (allowUpwind && !arg){
    int signizq,signder;
    double l1i=eigenvalueTab(1,pri[j],eos), l2i=eigenvalueTab(2,pri[j],eos), l3i=eigenvalueTab(3,pri[j],eos), l1d=eigenvalueTab(1,pri[j+1],eos), l2d=eigenvalueTab(2,pri[j+1],eos), l3d=eigenvalueTab(3,pri[j+1],eos);
    signizq=SIGN(l1i)+SIGN(l2i)+SIGN(l3i);
    signder=SIGN(l1d)+SIGN(l2d)+SIGN(l3d);
    if(l1i*l1d>0 && l2i*l2d>0 && l3i*l3d>0 && (signizq==3 || signizq==-3) && (signder==3 || signder==-3)){ //just if possible
      double *phiplustemp, *phimintemp;
      if((phiplustemp = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory phipt\n"); exit(41);}
      if((phimintemp = (double*)malloc(dim*sizeof(double))) == NULL) {printf("Error allocating memory phimt\n"); exit(41);}

      for(i=0;i<2*method+1;i++){
        for(k=0;k<dim;k++){
          f(valuef, u[j-method+i],pri[j-method+i]);
          m[i][k]=dot(l[k],valuef,dim);
        }
      }
      reconstruct(method,m,method,dim,dx,phiplustemp,empty);

      for(i=0;i<2*method+1;i++){
        for(k=0;k<dim;k++){
          f(valuef, u[j-method+i+1],pri[j-method+i+1]);
          m[i][k]=dot(l[dim+k],valuef,dim);
        }
      }
      reconstruct(method,m,method,dim,dx,empty,phimintemp);

      for (k=0;k<3;k++){
        if(eigenvalueTab(k+1,pri[j],eos)>0) {
          phiplus[k]=phiplustemp[k];
          phimin[k]=0.0;
        }
        else{
          phiplus[k]=0.0;
          phimin[k]=phimintemp[k];
        }
      }
      cleanM(m,(2*method+1)); free(valuef); free(phimintemp); free(phiplustemp); free(empty);
      return 0;
    }
  }
  
  double *alpha; double max1, max2, max3; int state, statem1;
  if((alpha = (double*)malloc((4+2*(method-1))*sizeof(double))) == NULL) {printf("Error allocating memory alpha\n"); exit(41);}
  for(i=0;i<4+2*(method-1);i++){
    state=-method+i; if (j+state>=N+5) statem1=state; else statem1=state+1;
     max1=fmax(fabs(eigenvalueTab(1,pri[j+state],eos)),fabs(eigenvalueTab(1,pri[j+statem1],eos)));
     max2=fmax(fabs(eigenvalueTab(2,pri[j+state],eos)),fabs(eigenvalueTab(2,pri[j+statem1],eos)));
     max3=fmax(fabs(eigenvalueTab(3,pri[j+state],eos)),fabs(eigenvalueTab(3,pri[j+statem1],eos)));
     max1=fmax(max1,max2);
     alpha[i]=fmax(max1,max3);
     if(arg) alpha[i]=sqrt(0.5*(1+alpha[i]*alpha[i]));
  }
  for(i=1;i<4+2*(method-1);i++) alpha[0]=fmax(alpha[0],alpha[i]);

  for(i=0;i<2*method+1;i++){
    for(k=0;k<dim;k++){
      f(valuef, u[j-method+i],pri[j-method+i]);
      m[i][k]=0.5*(dot(l[k],valuef,dim)+alpha[0]*dot(l[k],u[j-method+i],dim));
    }
  }
  reconstruct(method, m, method, dim, dx, phiplus, empty);

  for(i=0;i<2*method+1;i++){
    for(k=0;k<dim;k++){
      f(valuef, u[j-method+i+1],pri[j-method+i+1]);
      m[i][k]=0.5*(dot(l[dim+k],valuef,dim)-alpha[0]*dot(l[dim+k],u[j-method+i+1],dim));
    }
  }
  reconstruct(method, m, method, dim, dx, empty, phimin);

  cleanM(m, 2*method+1); free(alpha); free(empty); free(valuef);
  return 0;
  
}
/*-------------------------------------*/

void FM_PP(struct PP eos, int allowUpwind,int method,double dx,double* result,int j, double **u, double** pri, double**l,double**r,int N){
 double *phiplus,*phimin;
  if( (phimin=(double*)malloc(3*sizeof(double)))==NULL ||(phiplus=(double*)malloc(3*sizeof(double)))==NULL  ){ printf("Error allocating memory for psi\n"); exit(0); }

  int k;
  double *aux1,*aux2,*aux3;
  if( (aux1=(double*)malloc(3*sizeof(double)))==NULL ||(aux2=(double*)malloc(3*sizeof(double)))==NULL  ||(aux3=(double*)malloc(3*sizeof(double)))==NULL  ){
    printf("Error allocating auxiliar memory\n");
    exit(0);
  }
  //singular interface: apply to every k, as it's singular for (u,v)=>for k=1,3. else: (global) LF or upwind
  int nonl=0;
  if(nonlfactorsignPP(pri[j][0], eos)*nonlfactorsignPP(pri[j+1][0], eos)<=0) nonl=1; //calculated at 1st order
  scalarsPP(eos,allowUpwind,method,dx,phiplus,phimin,j,u,pri,l,nonl,N);

  for(k=0;k<3;k++){
    vectorbyn(aux1,phiplus[k],r[k],3);
    vectorbyn(aux2,phimin[k],r[k+3],3);
    sumvecs(aux3,aux1,aux2,3);
    sumvecs(result,result,aux3,3);
  }
  
  clean(aux1);clean(aux2);clean(aux3); free(phiplus); free(phimin);
}
void FM_TASPP(struct TASPP eos, int allowUpwind,int method,double dx,double* result,int j, double **u, double** pri, double**l,double**r,int N){
 double *phiplus,*phimin;
  if( (phimin=(double*)malloc(3*sizeof(double)))==NULL ||(phiplus=(double*)malloc(3*sizeof(double)))==NULL  ){ printf("Error allocating memory for psi\n"); exit(0); }

  int k;
  double *aux1,*aux2,*aux3;
  if( (aux1=(double*)malloc(3*sizeof(double)))==NULL ||(aux2=(double*)malloc(3*sizeof(double)))==NULL  ||(aux3=(double*)malloc(3*sizeof(double)))==NULL  ){
    printf("Error allocating auxiliar memory\n");
    exit(0);
  }
  //singular interface: apply to every k, as it's singular for (u,v)=>for k=1,3. else: (global) LF or upwind
  int nonl=0;
  if(nonlfactorsignTASPP(pri[j][0], eos)*nonlfactorsignTASPP(pri[j+1][0], eos)<=0) nonl=1; //calculated at 1st order
  scalarsTASPP(eos,allowUpwind,method,dx,phiplus,phimin,j,u,pri,l,nonl,N);

  for(k=0;k<3;k++){
    vectorbyn(aux1,phiplus[k],r[k],3);
    vectorbyn(aux2,phimin[k],r[k+3],3);
    sumvecs(aux3,aux1,aux2,3);
    sumvecs(result,result,aux3,3);
  }
  
  clean(aux1);clean(aux2);clean(aux3); free(phiplus); free(phimin);
}
void FM_Tab(struct Tab eos, int allowUpwind,int method,double dx,double* result,int j, double **u, double** pri, double**l,double**r,int N){
 double *phiplus,*phimin;
  if( (phimin=(double*)malloc(3*sizeof(double)))==NULL ||(phiplus=(double*)malloc(3*sizeof(double)))==NULL  ){ printf("Error allocating memory for psi\n"); exit(0); }

  int k;
  double *aux1,*aux2,*aux3;
  if( (aux1=(double*)malloc(3*sizeof(double)))==NULL ||(aux2=(double*)malloc(3*sizeof(double)))==NULL  ||(aux3=(double*)malloc(3*sizeof(double)))==NULL  ){
    printf("Error allocating auxiliar memory\n");
    exit(0);
  }
  //singular interface: apply to every k, as it's singular for (u,v)=>for k=1,3. else: (global) LF or upwind
  int nonl=0;
  if(nonlfactorsignTab(pri[j][0], eos)*nonlfactorsignTab(pri[j+1][0], eos)<=0) nonl=1; //calculated at 1st order
  scalarsTab(eos,allowUpwind,method,dx,phiplus,phimin,j,u,pri,l,nonl,N);

  for(k=0;k<3;k++){
    vectorbyn(aux1,phiplus[k],r[k],3);
    vectorbyn(aux2,phimin[k],r[k+3],3);
    sumvecs(aux3,aux1,aux2,3);
    sumvecs(result,result,aux3,3);
  }
  
  clean(aux1);clean(aux2);clean(aux3); free(phiplus); free(phimin);
}
/*----------------------------*/

void integratePP(struct PP eos, int allowUpwind, int method, int N, double**in, double **out, double **pri, double dx, double cfl, double**l, double **r){
  int dim=3; int nvarpri=3; //eran 4 por el pratio
  double *fluxL,*fluxR;
  if( (fluxL=(double*)malloc(dim*sizeof(double)))==NULL ||(fluxR=(double*)malloc(dim*sizeof(double)))==NULL){
    printf("Error allocating flux auxiliar memory\n");
    exit(0);
 }
  double *aux1,*aux2,*aux3;
  if( (aux1=(double*)malloc(dim*sizeof(double)))==NULL ||(aux2=(double*)malloc(dim*sizeof(double)))==NULL  ||(aux3=(double*)malloc(dim*sizeof(double)))==NULL  ){
    printf("Error allocating auxiliar memory\n");
    exit(0);
  }

  double **recons; int j;
  if ((recons=(double**) malloc(6*sizeof(double*)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}
  for (j=0;j<6;j++) if ((recons[j]=(double*) malloc(nvarpri*sizeof(double)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}

  double *left, *right;
  if((left=(double*)malloc(nvarpri*sizeof(double)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}
  if( (right=(double*)malloc(nvarpri*sizeof(double)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}
  
  for(j=3;j<N+3;j++){
    //reconstruction of interfaces
    if (j==3){ //first point
      memcpy(recons[0], pri[0], nvarpri*sizeof(double));
      memcpy(recons[1], pri[0], nvarpri*sizeof(double));
      reconstruct(method, pri, j, nvarpri, dx,right, left);
      memcpy(recons[2], left, nvarpri*sizeof(double));
      memcpy(recons[3], right, nvarpri*sizeof(double));
      reconstruct(method, pri, j+1, nvarpri, dx, right, left);
      memcpy(recons[4], left, nvarpri*sizeof(double));
      memcpy(recons[5], right, nvarpri*sizeof(double));
    }
    else if (j==N+2){ //last point
      displace(recons,nvarpri);
      memcpy(recons[4], pri[N+3], nvarpri*sizeof(double));
      memcpy(recons[5], pri[N+3], nvarpri*sizeof(double));
    }
    else{ //any other
      displace(recons,nvarpri);
      reconstruct(method, pri, j+1, nvarpri, dx, right, left);
      memcpy(recons[4], left, nvarpri*sizeof(double));
      memcpy(recons[5], right, nvarpri*sizeof(double));
    }
    
    //calculation of numerical fluxes
    memset(fluxL,0,dim*sizeof(double));memset(fluxR,0,dim*sizeof(double));
    eigenvectorsPP(eos,r,l,recons[3],recons[4]);
    FM_PP(eos, allowUpwind, method, dx, fluxL, j, in, pri, l, r, N);
    eigenvectorsPP(eos,r,l,recons[1],recons[2]);
    FM_PP(eos, allowUpwind, method, dx, fluxR, j-1, in, pri, l, r, N);

    //Apply fluxes
    vectorbyn(aux1,-1.0,fluxR,dim);
    sumvecs(aux2,fluxL,aux1,dim);
    vectorbyn(aux3,-cfl,aux2,dim);
    sumvecs(out[j],in[j],aux3,dim); 
  }
  
  //Free memory
  clean(fluxL);clean(fluxR);clean(aux1);clean(aux2);clean(aux3);
}
void integrateTASPP(struct TASPP eos, int allowUpwind, int method, int N, double**in, double **out, double **pri, double dx, double cfl, double**l, double **r){
  int dim=3; int nvarpri=3; //eran 4 por el pratio
  double *fluxL,*fluxR;
  if( (fluxL=(double*)malloc(dim*sizeof(double)))==NULL ||(fluxR=(double*)malloc(dim*sizeof(double)))==NULL){
    printf("Error allocating flux auxiliar memory\n");
    exit(0);
 }
  double *aux1,*aux2,*aux3;
  if( (aux1=(double*)malloc(dim*sizeof(double)))==NULL ||(aux2=(double*)malloc(dim*sizeof(double)))==NULL  ||(aux3=(double*)malloc(dim*sizeof(double)))==NULL  ){
    printf("Error allocating auxiliar memory\n");
    exit(0);
  }

  double **recons; int j;
  if ((recons=(double**) malloc(6*sizeof(double*)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}
  for (j=0;j<6;j++) if ((recons[j]=(double*) malloc(nvarpri*sizeof(double)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}

  double *left, *right;
  if((left=(double*)malloc(nvarpri*sizeof(double)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}
  if( (right=(double*)malloc(nvarpri*sizeof(double)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}
  
  for(j=3;j<N+3;j++){
    //reconstruction of interfaces
    if (j==3){ //first point
      memcpy(recons[0], pri[0], nvarpri*sizeof(double));
      memcpy(recons[1], pri[0], nvarpri*sizeof(double));
      reconstruct(method, pri, j, nvarpri, dx,right, left);
      memcpy(recons[2], left, nvarpri*sizeof(double));
      memcpy(recons[3], right, nvarpri*sizeof(double));
      reconstruct(method, pri, j+1, nvarpri, dx, right, left);
      memcpy(recons[4], left, nvarpri*sizeof(double));
      memcpy(recons[5], right, nvarpri*sizeof(double));
    }
    else if (j==N+2){ //last point
      displace(recons,nvarpri);
      memcpy(recons[4], pri[N+3], nvarpri*sizeof(double));
      memcpy(recons[5], pri[N+3], nvarpri*sizeof(double));
    }
    else{ //any other
      displace(recons,nvarpri);
      reconstruct(method, pri, j+1, nvarpri, dx, right, left);
      memcpy(recons[4], left, nvarpri*sizeof(double));
      memcpy(recons[5], right, nvarpri*sizeof(double));
    }
    
    //calculation of numerical fluxes
    memset(fluxL,0,dim*sizeof(double));memset(fluxR,0,dim*sizeof(double));
    eigenvectorsTASPP(eos,r,l,recons[3],recons[4]);
    FM_TASPP(eos, allowUpwind, method, dx, fluxL, j, in, pri, l, r, N);
    eigenvectorsTASPP(eos,r,l,recons[1],recons[2]);
    FM_TASPP(eos, allowUpwind, method, dx, fluxR, j-1, in, pri, l, r, N);

    //Apply fluxes
    vectorbyn(aux1,-1.0,fluxR,dim);
    sumvecs(aux2,fluxL,aux1,dim);
    vectorbyn(aux3,-cfl,aux2,dim);
    sumvecs(out[j],in[j],aux3,dim); 
  }
  
  //Free memory
  clean(fluxL);clean(fluxR);clean(aux1);clean(aux2);clean(aux3);
}
void integrateTab(struct Tab eos, int allowUpwind, int method, int N, double**in, double **out, double **pri, double dx, double cfl, double**l, double **r){
  int dim=3; int nvarpri=3; //eran 4 por el pratio
  double *fluxL,*fluxR;
  if( (fluxL=(double*)malloc(dim*sizeof(double)))==NULL ||(fluxR=(double*)malloc(dim*sizeof(double)))==NULL){
    printf("Error allocating flux auxiliar memory\n");
    exit(0);
 }
  double *aux1,*aux2,*aux3;
  if( (aux1=(double*)malloc(dim*sizeof(double)))==NULL ||(aux2=(double*)malloc(dim*sizeof(double)))==NULL  ||(aux3=(double*)malloc(dim*sizeof(double)))==NULL  ){
    printf("Error allocating auxiliar memory\n");
    exit(0);
  }

  double **recons; int j;
  if ((recons=(double**) malloc(6*sizeof(double*)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}
  for (j=0;j<6;j++) if ((recons[j]=(double*) malloc(nvarpri*sizeof(double)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}

  double *left, *right;
  if((left=(double*)malloc(nvarpri*sizeof(double)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}
  if( (right=(double*)malloc(nvarpri*sizeof(double)))==NULL) {printf("Error allocating reconstruction memory\n"); exit(5);}
  
  for(j=3;j<N+3;j++){
    //reconstruction of interfaces
    if (j==3){ //first point
      memcpy(recons[0], pri[0], nvarpri*sizeof(double));
      memcpy(recons[1], pri[0], nvarpri*sizeof(double));
      reconstruct(method, pri, j, nvarpri, dx,right, left);
      memcpy(recons[2], left, nvarpri*sizeof(double));
      memcpy(recons[3], right, nvarpri*sizeof(double));
      reconstruct(method, pri, j+1, nvarpri, dx, right, left);
      memcpy(recons[4], left, nvarpri*sizeof(double));
      memcpy(recons[5], right, nvarpri*sizeof(double));
    }
    else if (j==N+2){ //last point
      displace(recons,nvarpri);
      memcpy(recons[4], pri[N+3], nvarpri*sizeof(double));
      memcpy(recons[5], pri[N+3], nvarpri*sizeof(double));
    }
    else{ //any other
      displace(recons,nvarpri);
      reconstruct(method, pri, j+1, nvarpri, dx, right, left);
      memcpy(recons[4], left, nvarpri*sizeof(double));
      memcpy(recons[5], right, nvarpri*sizeof(double));
    }
    
    //calculation of numerical fluxes
    memset(fluxL,0,dim*sizeof(double));memset(fluxR,0,dim*sizeof(double));
    eigenvectorsTab(eos,r,l,recons[3],recons[4]);
    FM_Tab(eos, allowUpwind, method, dx, fluxL, j, in, pri, l, r, N);
    eigenvectorsTab(eos,r,l,recons[1],recons[2]);
    FM_Tab(eos, allowUpwind, method, dx, fluxR, j-1, in, pri, l, r, N);

    //Apply fluxes
    vectorbyn(aux1,-1.0,fluxR,dim);
    sumvecs(aux2,fluxL,aux1,dim);
    vectorbyn(aux3,-cfl,aux2,dim);
    sumvecs(out[j],in[j],aux3,dim); 
  }
  
  //Free memory
  clean(fluxL);clean(fluxR);clean(aux1);clean(aux2);clean(aux3);
}
/*----------------------------*/

void rk3(int usePP, int useTASPP, void *eosGeneric, double *x, int allowUpwind,int method, double tf, int N,double**in,double **out, double **pri,double dx, double dt,  double cfl,double**l,double **r){
  int dim=3;
  double **temp2;
  temp2=(double**) malloc((N+6)*sizeof(double*));
  int i; for (i=0; i<N+6; i++) {temp2[i]=(double*)malloc(dim*sizeof(double));}

  int Time=(int)floor(tf/dt);
  printf("Time iterations: %d\n", Time);
  int n,k;
  double **swap;

  if(usePP){
    printf("Integrating a PP EoS\n");
    struct PP eos=*((struct PP *)eosGeneric);
    for(n=0;n<=Time;n++){
     // if(!(n%10)){printf("n=%d, %e %e\n",n,pri[512][0],pri[512][2]);}
      integratePP(eos, allowUpwind,method,N,in,out,pri,dx,cfl,l,r); //in=u^n  out=u1=S(u^n)
      recoverPP(eos,out,pri,N); 
      
      //Neumann homogeneous boundary conditions
      for (i=0;i<3;i++) {memcpy(out[i],out[3], dim*sizeof(double)); memcpy(pri[i],pri[3], dim*sizeof(double));} 
      for (i=N+3;i<N+6;i++) {memcpy(out[i],out[N+2], dim*sizeof(double)); memcpy(pri[i],pri[N+2], dim*sizeof(double));} 
    
      integratePP(eos, allowUpwind,method,N,out,temp2,pri,dx,cfl,l,r); //temp2=S(u1)
      for(i=3;i<N+3; i++) {for(k=0;k<dim;k++) out[i][k]=0.75*in[i][k]+0.25*temp2[i][k]; }     
      recoverPP(eos,out,pri,N); 
      
      //Neumann homogeneous boundary conditions
      for (i=0;i<3;i++) {memcpy(out[i],out[3], dim*sizeof(double)); memcpy(pri[i],pri[3], dim*sizeof(double));} 
      for (i=N+3;i<N+6;i++) {memcpy(out[i],out[N+2], dim*sizeof(double)); memcpy(pri[i],pri[N+2], dim*sizeof(double));} 
    
      integratePP(eos, allowUpwind,method,N,out,temp2,pri,dx,cfl,l,r);
      for(i=3;i<N+3; i++) {for(k=0;k<dim;k++) out[i][k]=1.0/3*in[i][k]+2.0/3*temp2[i][k];} 
      recoverPP(eos,out,pri,N);
      
      //Neumann homogeneous boundary conditions
      for (i=0;i<3;i++) {memcpy(out[i],out[3], dim*sizeof(double)); memcpy(pri[i],pri[3], dim*sizeof(double));} 
      for (i=N+3;i<N+6;i++) {memcpy(out[i],out[N+2], dim*sizeof(double)); memcpy(pri[i],pri[N+2], dim*sizeof(double));} 

      swap= in; in=out; out=swap;
    }
  }

  else if(useTASPP){
    printf("Integrating a TASPP EoS\n");

    struct TASPP eos=*((struct TASPP *)eosGeneric);
    for(n=0;n<=Time;n++){
      integrateTASPP(eos, allowUpwind,method,N,in,out,pri,dx,cfl,l,r); //in=u^n  out=u1=S(u^n)
      recoverTASPP(eos,out,pri,N); 
      
      //Neumann homogeneous boundary conditions
      for (i=0;i<3;i++) {memcpy(out[i],out[3], dim*sizeof(double)); memcpy(pri[i],pri[3], dim*sizeof(double));} 
      for (i=N+3;i<N+6;i++) {memcpy(out[i],out[N+2], dim*sizeof(double)); memcpy(pri[i],pri[N+2], dim*sizeof(double));} 
    
      integrateTASPP(eos, allowUpwind,method,N,out,temp2,pri,dx,cfl,l,r); //temp2=S(u1)
      for(i=3;i<N+3; i++) {for(k=0;k<dim;k++) out[i][k]=0.75*in[i][k]+0.25*temp2[i][k]; }     
      recoverTASPP(eos,out,pri,N); 
      
      //Neumann homogeneous boundary conditions
      for (i=0;i<3;i++) {memcpy(out[i],out[3], dim*sizeof(double)); memcpy(pri[i],pri[3], dim*sizeof(double));} 
      for (i=N+3;i<N+6;i++) {memcpy(out[i],out[N+2], dim*sizeof(double)); memcpy(pri[i],pri[N+2], dim*sizeof(double));} 
    
      integrateTASPP(eos, allowUpwind,method,N,out,temp2,pri,dx,cfl,l,r);
      for(i=3;i<N+3; i++) {for(k=0;k<dim;k++) out[i][k]=1.0/3*in[i][k]+2.0/3*temp2[i][k];} 
      recoverTASPP(eos,out,pri,N);
      
      //Neumann homogeneous boundary conditions
      for (i=0;i<3;i++) {memcpy(out[i],out[3], dim*sizeof(double)); memcpy(pri[i],pri[3], dim*sizeof(double));} 
      for (i=N+3;i<N+6;i++) {memcpy(out[i],out[N+2], dim*sizeof(double)); memcpy(pri[i],pri[N+2], dim*sizeof(double));} 
      swap= in; in=out; out=swap;
    }
  }
  else{  //Tab guaranteed from main
    printf("Integrating a Tabulated EoS\n");
    struct Tab eos=*((struct Tab *)eosGeneric);
    for(n=0;n<=Time;n++){
      integrateTab(eos, allowUpwind,method,N,in,out,pri,dx,cfl,l,r); //in=u^n  out=u1=S(u^n)
      recoverTab(eos,out,pri,N); 
      
      //Neumann homogeneous boundary conditions
      for (i=0;i<3;i++) {memcpy(out[i],out[3], dim*sizeof(double)); memcpy(pri[i],pri[3], dim*sizeof(double));} 
      for (i=N+3;i<N+6;i++) {memcpy(out[i],out[N+2], dim*sizeof(double)); memcpy(pri[i],pri[N+2], dim*sizeof(double));} 
    
      integrateTab(eos, allowUpwind,method,N,out,temp2,pri,dx,cfl,l,r); //temp2=S(u1)
      for(i=3;i<N+3; i++) {for(k=0;k<dim;k++) out[i][k]=0.75*in[i][k]+0.25*temp2[i][k]; }     
      recoverTab(eos,out,pri,N); 
      
      //Neumann homogeneous boundary conditions
      for (i=0;i<3;i++) {memcpy(out[i],out[3], dim*sizeof(double)); memcpy(pri[i],pri[3], dim*sizeof(double));} 
      for (i=N+3;i<N+6;i++) {memcpy(out[i],out[N+2], dim*sizeof(double)); memcpy(pri[i],pri[N+2], dim*sizeof(double));} 
    
      integrateTab(eos, allowUpwind,method,N,out,temp2,pri,dx,cfl,l,r);
      for(i=3;i<N+3; i++) {for(k=0;k<dim;k++) out[i][k]=1.0/3*in[i][k]+2.0/3*temp2[i][k];} 
      recoverTab(eos,out,pri,N);
      
      //Neumann homogeneous boundary conditions
      for (i=0;i<3;i++) {memcpy(out[i],out[3], dim*sizeof(double)); memcpy(pri[i],pri[3], dim*sizeof(double));} 
      for (i=N+3;i<N+6;i++) {memcpy(out[i],out[N+2], dim*sizeof(double)); memcpy(pri[i],pri[N+2], dim*sizeof(double));} 
      
      swap= in; in=out; out=swap;
    }
  }
  cleanM(temp2,N+6);
}



