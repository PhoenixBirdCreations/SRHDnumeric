#include "HLLrel.h"
#include "read_par_eos.h"
#include "read_par_ic.h"


int main(){
  int usePP=0, useTab=0, useTASPP=0;

  struct PP eosPP; struct TASPP eosTASPP; struct Tab eosTab;
  read_par_eos(&eosPP, &eosTASPP, &eosTab, "eos.par", &usePP, &useTab, &useTASPP);

  double *enteringInfo;
  if((enteringInfo=(double*)malloc(8*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
  char *filename;
  int *intInfo;
  if((intInfo=(int*)malloc(3*sizeof(int)))==NULL){printf("error reserving memory\n"); exit(1);}

  read_par_ic("ic.par", &filename, enteringInfo, intInfo);
  
  double rho_L, v_L, rho_R, v_R; double a, b; double cfl; double tf; 
  int N; int allowUpwind; int method;

  rho_L=enteringInfo[0]; v_L=enteringInfo[1]; rho_R=enteringInfo[2]; v_R=enteringInfo[3];
  a=enteringInfo[4]; b=enteringInfo[5];
  N=intInfo[0];
  cfl=enteringInfo[6]; tf=enteringInfo[7];
  allowUpwind=intInfo[1]; method=intInfo[2];
  printf("Output file: %s\n",filename);
  free(enteringInfo); free(intInfo);


  //------ End of reading -------//
  double dt=cfl*(b-a)/(N-1), dx=(b-a)/(N-1);
  int i; double Pl,Pr,hl,hr;

  //ensure selection of an EoS type
  if((useTASPP+usePP+useTab)!=1){printf("Select one and only one type of EoS\n"); exit(1);}

  //define the existence of eigenvectors
  double **r; double **l;
  if((r=(double**)malloc(6*sizeof(double*)))==NULL){
      printf("failed to allocate memory r gordo\n");return -1;}
  if((l=(double**)malloc(6*sizeof(double*)))==NULL){
      printf("failed to allocate memory l gordo\n");return -1;}
  for (i=0;i<6;i++){
      if((r[i]=(double*)malloc(3*sizeof(double)))==NULL){
      printf("failed to allocate memory r_comp\n");return -1;}
      if((l[i]=(double*)malloc(3*sizeof(double)))==NULL){
      printf("failed to allocate memory l_comp\n");return -1;}
  }  
    
  
  //define the existance of x,pri,u(matrix),temp(matrix)
  double *x; double **u; double **temp; double **pri;
  if((x=(double*)malloc(N*sizeof(double)))==NULL){
      printf("failed to allocate memory x\n");return -1;}
  if((u=(double**)malloc((N+6)*sizeof(double*)))==NULL){
      printf("failed to allocate memory u\n");return -1;}  
  if((pri=(double**)malloc((N+6)*sizeof(double*)))==NULL){
      printf("failed to allocate memory pri\n");return -1;} 
  if((temp=(double**)malloc((N+6)*sizeof(double*)))==NULL){
      printf("failed to allocate memory temp\n");return -1;}
  for(i=0;i<N+6;i++){
    if((u[i]=(double*)malloc(3*sizeof(double)))==NULL){
      printf("failed to allocate memory u_comp\n");return -1;}
    if((pri[i]=(double*)malloc(3*sizeof(double)))==NULL){
      printf("failed to allocate memory pri_comp\n");return -1;}
    if((temp[i]=(double*)malloc(3*sizeof(double)))==NULL){
      printf("failed to allocate memory temp_comp\n");return -1;}
  }

  //initialize and integrate using the selected type of EoS
  if (usePP){ //EoS is piecewise polytropic
    printf("Initializing\n");
    ciPP(eosPP,v_L,rho_L,v_R,rho_R,&Pl,&Pr,&hl,&hr); 
    initialize(x,u,pri,a,b,N,rho_L,v_L,Pl,hl,rho_R,v_R,Pr,hr);
    memcpy(temp[0],u[0],3*sizeof(double)); memcpy(temp[N+5],u[N+5],3*sizeof(double));
    memcpy(temp[1],u[1],3*sizeof(double)); memcpy(temp[N+4],u[N+4],3*sizeof(double));
    memcpy(temp[2],u[2],3*sizeof(double)); memcpy(temp[N+3],u[N+3],3*sizeof(double));
    
    printf("integrating\n");
    time_t start=clock();
    rk3(usePP,useTASPP,&eosPP,x,allowUpwind,method,tf,N,u,temp,pri,dx,dt,cfl,l,r);
    double total_time=(double) (clock()-start)/CLOCKS_PER_SEC;
    printf("Time passed: %fs, equivalent to %fm or to %fh\n", total_time, total_time/60, total_time/3600);
    
    if(((int)(floor(tf/dt))+1)%2==0) {double **swap; swap=u; u=temp; temp=swap;}
    writePP(eosPP,x,u,pri,N,filename);
    cleanPP(&eosPP);
  }

  else if (useTASPP){//EoS is TASPP
    printf("Initializing\n");
    ciTASPP(eosTASPP,v_L,rho_L,v_R,rho_R,&Pl,&Pr,&hl,&hr);
    initialize(x,u,pri,a,b,N,rho_L,v_L,Pl,hl,rho_R,v_R,Pr,hr);
    memcpy(temp[0],u[0],3*sizeof(double)); memcpy(temp[N+5],u[N+5],3*sizeof(double));
    memcpy(temp[1],u[1],3*sizeof(double)); memcpy(temp[N+4],u[N+4],3*sizeof(double));
    memcpy(temp[2],u[2],3*sizeof(double)); memcpy(temp[N+3],u[N+3],3*sizeof(double));
    
    printf("integrating\n");
    time_t start=clock();
    rk3(usePP,useTASPP,&eosTASPP,x,allowUpwind,method,tf,N,u,temp,pri,dx,dt,cfl,l,r);
    double total_time=(double) (clock()-start)/CLOCKS_PER_SEC;
    printf("Time passed: %fs, equivalent to %fm or to %fh\n", total_time, total_time/60, total_time/3600);
    
    if(((int)(floor(tf/dt))+1)%2==0) {double **swap; swap=u; u=temp; temp=swap;}
    writeTASPP(eosTASPP,x,u,pri,N,filename);
    cleanTASPP(&eosTASPP);
  }


  else if (useTab){//EoS is tabulated
    printf("Initializing\n");
    ciTAB(eosTab,v_L,rho_L,v_R,rho_R,&Pl,&Pr,&hl,&hr);
    initialize(x,u,pri,a,b,N,rho_L,v_L,Pl,hl,rho_R,v_R,Pr,hr);
    memcpy(temp[0],u[0],3*sizeof(double)); memcpy(temp[N+5],u[N+5],3*sizeof(double));
    memcpy(temp[1],u[1],3*sizeof(double)); memcpy(temp[N+4],u[N+4],3*sizeof(double));
    memcpy(temp[2],u[2],3*sizeof(double)); memcpy(temp[N+3],u[N+3],3*sizeof(double));
      
    printf("integrating\n");
    time_t start=clock();
    rk3(usePP,useTASPP,&eosTab,x,allowUpwind,method,tf,N,u,temp,pri,dx,dt,cfl,l,r);
    double total_time=(double) (clock()-start)/CLOCKS_PER_SEC;
    printf("Time passed: %fs, equivalent to %fm or to %fh\n", total_time, total_time/60, total_time/3600);
    
    if(((int)(floor(tf/dt))+1)%2==0) {double **swap; swap=u; u=temp; temp=swap;}
    writeTab(eosTab,x,u,pri,N,filename);
    cleanTab(&eosTab); 
  }

  //free memory and exit
  clean(x);cleanM(r,6);cleanM(l,6);cleanM(u,N+6);cleanM(temp,N+6); cleanM(pri,N+6);
  return 0;

}

