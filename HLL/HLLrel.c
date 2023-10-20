#include "HLLrel.h"

#define SIGN(x) (((x)>0) ? 1 : -1)

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN3(a,b,c) MIN(MIN(a,b),c)
#define MIN4(a,b,c,d) MIN(MIN(a,b),MIN(c,d))
#define MIN5(a,b,c,d,e) MIN(MIN4(a,b,c,d),e)
#define MIN6(a,b,c,d,e,f) MIN(MIN4(a,b,c,d),MIN(e,f))
#define MIN7(a,b,c,d,e,f,g) MIN(MIN4(a,b,c,d),MIN3(e,f,g))
#define MIN8(a,b,c,d,e,f,g,h) MIN(MIN4(a,b,c,d),MIN4(e,f,g,h))
#define MAX3(a,b,c) MAX(MAX(a,b),c)
#define MAX4(a,b,c,d) MAX(MAX(a,b),MAX(c,d))
#define MAX5(a,b,c,d,e) MAX(MAX4(a,b,c,d),e)
#define MAX6(a,b,c,d,e,f) MAX(MAX4(a,b,c,d),MAX(e,f))
#define MAX7(a,b,c,d,e,f,g) MAX(MAX4(a,b,c,d),MAX3(e,f,g))
#define MAX8(a,b,c,d,e,f,g,h) MAX(MAX4(a,b,c,d),MAX4(e,f,g,h))
double max7(double a,double b,double c,double d,double e,double f, double g)
{
  return MAX8(0., a,b,c,d,e,f,g);
}
double min7(double a,double b,double c,double d,double e,double f, double g)
{
  return MIN8(0., a,b,c,d,e,f,g);
}

void FHLL_PP(struct PP eos, double* result,int j, double **u, double** pri, double *reconL, double *reconR){

  double max1,min1;
         
  double *ful,*fur;
  if( (ful=(double*)malloc(3*sizeof(double)))==NULL ||(fur=(double*)malloc(3*sizeof(double)))==NULL){
    printf("Error allocating auxiliar memory\n");
    exit(0);
  }

  f(ful, u[j], reconL); 
  f(fur, u[j+1], reconR); 
  
  double l1i, l1d, l2d,l2i,l3i,l3d;
  l1i=eigenvaluePP(1, pri[j], eos);l1d=eigenvaluePP(1,pri[j+1], eos);
  l2i=eigenvaluePP(2,pri[j], eos);l2d=eigenvaluePP(2,pri[j+1], eos);
  l3i=eigenvaluePP(3,pri[j], eos);l3d=eigenvaluePP(3,pri[j+1], eos);
  max1=MAX7(0.0, l1i,l1d,l2i,l2d,l3i,l3d);
  min1=MIN7(0.0, l1i,l1d,l2i,l2d,l3i,l3d);
  //like this is upwind interface, with speed estimated

  //sonic interface
  if (max1<=0.0){
    max1=MIN(0.5, fabs(min1)); 
  }
  else if (min1>=0.0){
    min1=-MIN(0.5, max1);
  }

  //nonconvex interface
  if(nonlfactorsignPP(pri[j][0], eos)*nonlfactorsignPP(pri[j+1][0], eos)<0){
    (max1)=sqrt((1+(max1))*0.5);
    (min1)=-sqrt((1+fabs((min1)))*0.5);
  } 

  int k; double dem=1.0/(max1-min1);

  //construct high order converved quantities from high order primitives
  double DL, DR, SL, SR, tauL, tauR;
  double w, he, v, rho;
  v=reconL[1]; rho=reconL[0];
  w=1.0/sqrt(1-v*v);
  he=hPP(rho, eos);
  DL=rho*w;
  SL=he*rho*w*w*v;
  tauL=he*rho*w*w-reconL[2]-DL;
  v=reconR[1]; rho=reconR[0];
  w=1.0/sqrt(1-v*v);
  he=hPP(rho, eos);
  DR=rho*w;
  SR=he*rho*w*w*v;
  tauR=he*rho*w*w-reconR[2]-DR;

  //HLL flux
  k=0;
  result[k]=(max1*ful[k]-min1*fur[k]+max1*min1*(DR-DL))*dem;
  k=1;
  result[k]=(max1*ful[k]-min1*fur[k]+max1*min1*(SR-SL))*dem;
  k=2;
  result[k]=(max1*ful[k]-min1*fur[k]+max1*min1*(tauR-tauL))*dem;

  if(isnan(result[0])) {
    printf("nan in flux: %f %f %f %f %f %f\n",ful[0],ful[1],ful[2],fur[0],fur[1],fur[2]);
    printf("reconR %f %f %f, prim %f %f %f\n", reconR[0], reconR[1], reconR[2],pri[j+1][0],pri[j+1][1],pri[j+1][2]);
    exit(4);
  }
    
  clean(fur);clean(ful);
}
void FHLL_TASPP(struct TASPP eos, double* result,int j, double **u, double** pri, double *reconL, double *reconR){

  double max1,min1;
         
  double *ful,*fur;
  if( (ful=(double*)malloc(3*sizeof(double)))==NULL ||(fur=(double*)malloc(3*sizeof(double)))==NULL){
    printf("Error allocating auxiliar memory\n");
    exit(0);
  }

  f(ful, u[j], reconL); 
  f(fur, u[j+1], reconR); 
  
  double l1i, l1d, l2d,l2i,l3i,l3d;
  l1i=eigenvalueTASPP(1, pri[j], eos);l1d=eigenvalueTASPP(1,pri[j+1], eos);
  l2i=eigenvalueTASPP(2,pri[j], eos);l2d=eigenvalueTASPP(2,pri[j+1], eos);
  l3i=eigenvalueTASPP(3,pri[j], eos);l3d=eigenvalueTASPP(3,pri[j+1], eos);
  max1=MAX7(0.0, l1i,l1d,l2i,l2d,l3i,l3d);
  min1=MIN7(0.0, l1i,l1d,l2i,l2d,l3i,l3d);
  //like this is upwind interface, with speed estimated

  //sonic interface
  if (max1<=0.0){
    max1=MIN(0.5, fabs(min1)); 
  }
  else if (min1>=0.0){
    min1=-MIN(0.5, max1);
  }

  //nonconvex interface
  if(nonlfactorsignTASPP(pri[j][0], eos)*nonlfactorsignTASPP(pri[j+1][0], eos)<0){
    (max1)=sqrt((1+(max1))*0.5);
    (min1)=-sqrt((1+fabs((min1)))*0.5);
  } 

  int k; double dem=1.0/(max1-min1);

  //construct high order converved quantities from high order primitives
  double DL, DR, SL, SR, tauL, tauR;
  double w, he, v, rho;
  v=reconL[1]; rho=reconL[0];
  w=1.0/sqrt(1-v*v);
  he=hTASPP(rho, eos);
  DL=rho*w;
  SL=he*rho*w*w*v;
  tauL=he*rho*w*w-reconL[2]-DL;
  v=reconR[1]; rho=reconR[0];
  w=1.0/sqrt(1-v*v);
  he=hTASPP(rho, eos);
  DR=rho*w;
  SR=he*rho*w*w*v;
  tauR=he*rho*w*w-reconR[2]-DR;

  //HLL flux
  k=0;
  result[k]=(max1*ful[k]-min1*fur[k]+max1*min1*(DR-DL))*dem;
  k=1;
  result[k]=(max1*ful[k]-min1*fur[k]+max1*min1*(SR-SL))*dem;
  k=2;
  result[k]=(max1*ful[k]-min1*fur[k]+max1*min1*(tauR-tauL))*dem;

  if(isnan(result[0])) {
    printf("nan in flux: %f %f %f %f %f %f\n",ful[0],ful[1],ful[2],fur[0],fur[1],fur[2]);
    printf("reconR %f %f %f, prim %f %f %f\n", reconR[0], reconR[1], reconR[2],pri[j+1][0],pri[j+1][1],pri[j+1][2]);
    exit(4);
  }
    
  clean(fur);clean(ful);
}
void FHLL_TAB(struct Tab eos, double* result,int j, double **u, double** pri, double *reconL, double *reconR){

  double max1,min1;
         
  double *ful,*fur;
  if( (ful=(double*)malloc(3*sizeof(double)))==NULL ||(fur=(double*)malloc(3*sizeof(double)))==NULL){
    printf("Error allocating auxiliar memory\n");
    exit(0);
  }

  f(ful, u[j], reconL); 
  f(fur, u[j+1], reconR); 
  
  double l1i, l1d, l2d,l2i,l3i,l3d;
  l1i=eigenvalueTab(1, pri[j], eos);l1d=eigenvalueTab(1,pri[j+1], eos);
  l2i=eigenvalueTab(2,pri[j], eos);l2d=eigenvalueTab(2,pri[j+1], eos);
  l3i=eigenvalueTab(3,pri[j], eos);l3d=eigenvalueTab(3,pri[j+1], eos);
  max1=MAX7(0.0, l1i,l1d,l2i,l2d,l3i,l3d);
  min1=MIN7(0.0, l1i,l1d,l2i,l2d,l3i,l3d);
  //like this is upwind interface, with speed estimated

  //sonic interface
  if (max1<=0.0){
    max1=MIN(0.5, fabs(min1)); 
  }
  else if (min1>=0.0){
    min1=-MIN(0.5, max1);
  }

  //nonconvex interface
  if(nonlfactorsignTab(pri[j][0], eos)*nonlfactorsignTab(pri[j+1][0], eos)<0){
    (max1)=sqrt((1+(max1))*0.5);
    (min1)=-sqrt((1+fabs((min1)))*0.5);
  } 

  int k; double dem=1.0/(max1-min1);

  //construct high order converved quantities from high order primitives
  double DL, DR, SL, SR, tauL, tauR;
  double w, he, v, rho;
  v=reconL[1]; rho=reconL[0];
  w=1.0/sqrt(1-v*v);
  he=hTab(rho, eos);
  DL=rho*w;
  SL=he*rho*w*w*v;
  tauL=he*rho*w*w-reconL[2]-DL;
  v=reconR[1]; rho=reconR[0];
  w=1.0/sqrt(1-v*v);
  he=hTab(rho, eos);
  DR=rho*w;
  SR=he*rho*w*w*v;
  tauR=he*rho*w*w-reconR[2]-DR;

  //HLL flux
  k=0;
  result[k]=(max1*ful[k]-min1*fur[k]+max1*min1*(DR-DL))*dem;
  k=1;
  result[k]=(max1*ful[k]-min1*fur[k]+max1*min1*(SR-SL))*dem;
  k=2;
  result[k]=(max1*ful[k]-min1*fur[k]+max1*min1*(tauR-tauL))*dem;

  if(isnan(result[0])) {
    printf("nan in flux: %f %f %f %f %f %f\n",ful[0],ful[1],ful[2],fur[0],fur[1],fur[2]);
    printf("reconR %f %f %f, prim %f %f %f\n", reconR[0], reconR[1], reconR[2],pri[j+1][0],pri[j+1][1],pri[j+1][2]);
    exit(4);
  }
    
  clean(fur);clean(ful);
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
    FHLL_PP(eos, fluxR, j, in, pri, recons[3],recons[4]);
    FHLL_PP(eos, fluxL, j-1, in, pri, recons[1],recons[2]);

    //Apply fluxes
    vectorbyn(aux1,-1.0,fluxL,dim);
    sumvecs(aux2,fluxR,aux1,dim);
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
    FHLL_TASPP(eos, fluxR, j, in, pri, recons[3],recons[4]);
    FHLL_TASPP(eos, fluxL, j-1, in, pri, recons[1],recons[2]);

    //Apply fluxes
    vectorbyn(aux1,-1.0,fluxL,dim);
    sumvecs(aux2,fluxR,aux1,dim);
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
    FHLL_TAB(eos, fluxR, j, in, pri, recons[3],recons[4]);
    FHLL_TAB(eos, fluxL, j-1, in, pri, recons[1],recons[2]);

    //Apply fluxes
    vectorbyn(aux1,-1.0,fluxL,dim);
    sumvecs(aux2,fluxR,aux1,dim);
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



