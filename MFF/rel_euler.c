#include "rel_euler.h"

double W(double x){
 return 1.0/sqrt(1-x*x); 
}

void f(double *re,double*u, double *pri){
 re[0]=u[0]*pri[1];re[1]=u[1]*pri[1]+pri[2];re[2]=u[1]-u[0]*pri[1];
}

//eigenvalues: 1_ (v-cs)/(1-vcs)  2_ v  3_ (v+cs)/(1+vcs)
double eigenvaluePP(int k, double *pri, struct PP eos){
 if(k==2)
   return pri[1];
 else{
   int s=k-2; 
   double v=pri[1], cs=sqrt(csPP(pri[0],eos));
   return (v+s*cs)/(1+s*v*cs);
 }
}
double eigenvalueTASPP(int k, double *pri, struct TASPP eos){
 if(k==2)
   return pri[1];
 else{
   int s=k-2; 
   double v=pri[1], cs=sqrt(csTASPP(pri[0],eos));
   return (v+s*cs)/(1+s*v*cs);
 }
}
double eigenvalueTab(int k, double *pri, struct Tab eos){
 if(k==2)
   return pri[1];
 else{
   int s=k-2; 
   double v=pri[1], cs=sqrt(csTab(pri[0],eos));
   return (v+s*cs)/(1+s*v*cs);
 }
}


void eigenvectorsPP(struct PP eos, double **r, double **l, double*pl, double*pr){ 
  int k,index;
  //left side
  double rho=pl[0],v=pl[1];
  double css=csPP(rho, eos);
  double c=sqrt(css);
  double he=hPP(rho,eos);
  double w=W(v);
  double kfancy=kfancyPP(rho,he,eos);
  //det
  double det=he*w*(kfancy-1)*(1-v*v)*2*c;
  double factor=1.0/det;

  for (k=1;k<4;k++){
   index=k-1;
   l[index][0]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v+(k-2)*he*w*c*(1-v*v));  //+-h^2/det*(k(v+-cs)-v-+hwcs(1-v^2))
   l[index][1]=(2-k)*factor*(1-kfancy*(1+(2-k)*v*c));  //+-h^2/det*(1-k*(1-+vcs))
   l[index][2]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v);  //+-h^2/det*(k(v-+cs)-v )
   
   r[index][0]=1.0;
   r[index][1]=he*w*(v+(k-2)*c);  //hw(v-+cs)
   r[index][2]=he*w*(1+(k-2)*v*c)-1;   //hw(1-+vcs)-1
   k++;
  }
  factor=w/(kfancy-1);
  l[1][0]=factor*(he-w); //w/(k-1)*(h-w)
  l[1][1]=w*v*factor; //w/(k-1)*wv
  l[1][2]=-w*factor;  //w/(k-1)*(-w)
  r[1][0]=kfancy/(he*w);//K/hw
  r[1][1]=v;  //v
  r[1][2]=1-r[1][0];  //1-k/hw

  
  //right side
  rho=pr[0],v=pr[1];
  css=csPP(rho, eos);
  c=sqrt(css);
  he=hPP(rho,eos);
  w=W(v);
  kfancy=kfancyPP(rho,he,eos);
  //det
  det=he*w*(kfancy-1)*(1-v*v)*2*c;
  factor=1.0/det;
  for (k=1;k<4;k++){
   index=3+k-1;
   l[index][0]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v+(k-2)*he*w*c*(1-v*v));
   l[index][1]=(2-k)*factor*(1-kfancy*(1+(2-k)*v*c));
   l[index][2]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v);
   r[index][0]=1.0;
   r[index][1]=he*w*(v+(k-2)*c);
   r[index][2]=he*w*(1+(k-2)*v*c)-1;
   k++;
  }
  factor=w/(kfancy-1);
  l[4][0]=factor*(he-w); 
  l[4][1]=w*v*factor; 
  l[4][2]=-w*factor;
  r[4][0]=kfancy/(he*w); 
  r[4][1]=v;
  r[4][2]=1-r[4][0]; 
}
void eigenvectorsTASPP(struct TASPP eos, double **r, double **l, double*pl, double*pr){ 
  int k,index;
  //left side
  double rho=pl[0],v=pl[1];
  double css=csTASPP(rho, eos);
  double c=sqrt(css);
  double he=hTASPP(rho,eos);
  double w=W(v);
  double kfancy=kfancyTASPP(rho,eos);
  //det
  double det=he*w*(kfancy-1)*(1-v*v)*2*c;
  double factor=1.0/det;

  for (k=1;k<4;k++){
   index=k-1;
   l[index][0]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v+(k-2)*he*w*c*(1-v*v));  //+-h^2/det*(k(v-+cs)-v+-hwcs(1-v^2))
   l[index][1]=(2-k)*factor*(1-kfancy*(1+(2-k)*v*c));  //+-h^2/det*(1-k*(1-+vcs))
   l[index][2]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v);  //+-h^2/det*(k(v-+cs)-v )
   
   r[index][0]=1.0;
   r[index][1]=he*w*(v+(k-2)*c);  //hw(v-+cs)
   r[index][2]=he*w*(1+(k-2)*v*c)-1;   //hw(1-+vcs)-1
   k++;
  }
  factor=w/(kfancy-1);
  l[1][0]=factor*(he-w); //w/(k-1)*(h-w)
  l[1][1]=w*v*factor; //w/(k-1)*wv
  l[1][2]=-w*factor;  //w/(k-1)*(-w)
  r[1][0]=kfancy/(he*w);//K/hw
  r[1][1]=v;  //v
  r[1][2]=1-r[1][0];  //1-k/hw

  
  //right side
  rho=pr[0],v=pr[1];
  css=csTASPP(rho, eos);
  c=sqrt(css);
  he=hTASPP(rho,eos);
  w=W(v);
  kfancy=kfancyTASPP(rho,eos);
  //det
  det=he*w*(kfancy-1)*(1-v*v)*2*c;
  factor=1.0/det;
  for (k=1;k<4;k++){
   index=3+k-1;
   l[index][0]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v+(k-2)*he*w*c*(1-v*v));
   l[index][1]=(2-k)*factor*(1-kfancy*(1+(2-k)*v*c));
   l[index][2]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v);
   r[index][0]=1.0;
   r[index][1]=he*w*(v+(k-2)*c);
   r[index][2]=he*w*(1+(k-2)*v*c)-1;
   k++;
  }
  factor=w/(kfancy-1);
  l[4][0]=factor*(he-w); 
  l[4][1]=w*v*factor; 
  l[4][2]=-w*factor;
  r[4][0]=kfancy/(he*w); 
  r[4][1]=v;
  r[4][2]=1-r[4][0]; 
}
void eigenvectorsTab(struct Tab eos, double **r, double **l, double*pl, double*pr){ 
  int k,index;
  //left side
  double rho=pl[0],v=pl[1];
  double css=csTab(rho, eos);
  double c=sqrt(css);
  double he=hTab(rho,eos);
  double w=W(v);
  double kfancy=kfancyTab(rho,eos);
  //det
  double det=he*w*(kfancy-1)*(1-v*v)*2*c;
  double factor=1.0/det;

  for (k=1;k<4;k++){
   index=k-1;
   l[index][0]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v+(k-2)*he*w*c*(1-v*v));//+-h^2/det*(k(v-+cs)-v+-hwcs(1-v^2))
   l[index][1]=(2-k)*factor*(1-kfancy*(1+(2-k)*v*c));  //+-h^2/det*(1-k*(1-+vcs))
   l[index][2]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v);  //+-h^2/det*(k(v-+cs)-v )
   
   r[index][0]=1.0;
   r[index][1]=he*w*(v+(k-2)*c);  //hw(v-+cs)
   r[index][2]=he*w*(1+(k-2)*v*c)-1;   //hw(1-+vcs)-1
   k++;
  }
  factor=w/(kfancy-1);
  l[1][0]=factor*(he-w); //w/(k-1)*(h-w)
  l[1][1]=w*v*factor; //w/(k-1)*wv
  l[1][2]=-w*factor;  //w/(k-1)*(-w)
  r[1][0]=kfancy/(he*w);//K/hw
  r[1][1]=v;  //v
  r[1][2]=1-r[1][0];  //1-k/hw

  
  //right side
  rho=pr[0],v=pr[1];
  css=csTab(rho, eos);
  c=sqrt(css);
  he=hTab(rho,eos);
  w=W(v);
  kfancy=kfancyTab(rho,eos);
  //det
  det=he*w*(kfancy-1)*(1-v*v)*2*c;
  factor=1.0/det;
  for (k=1;k<4;k++){
   index=3+k-1;
   l[index][0]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v+(k-2)*he*w*c*(1-v*v));
   l[index][1]=(2-k)*factor*(1-kfancy*(1+(2-k)*v*c));
   l[index][2]=(2-k)*factor*(kfancy*(v+(2-k)*c)-v);
   r[index][0]=1.0;
   r[index][1]=he*w*(v+(k-2)*c);
   r[index][2]=he*w*(1+(k-2)*v*c)-1;
   k++;
  }
  factor=w/(kfancy-1);
  l[4][0]=factor*(he-w); 
  l[4][1]=w*v*factor; 
  l[4][2]=-w*factor;
  r[4][0]=kfancy/(he*w); 
  r[4][1]=v;
  r[4][2]=1-r[4][0]; 
}


