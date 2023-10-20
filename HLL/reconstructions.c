#include "reconstructions.h"


void phm(double **result, double  **u, int j, int dim, double h){
    double dm, dp, d, alpha, tol, adm, adp, trans, nua, numa;
    int k;
    tol=h*h;

    for (k=0;k<dim; k++){
        dp=(u[j+1][k]-u[j][k])/h;
        dm=(u[j][k]-u[j-1][k])/h;
        adm=fabs(dm);
        adp=fabs(dp);
        trans=dp*dm;

        //selection
        if ((adm<=tol && adp<=tol) || trans<=0){
            d=0; alpha=0;
        }
        else{
            if (adm<=tol){
                d=2*dp*tol/(1+tol); alpha=2*(sqrt(2/(1+tol))-1);
            }
            else{
                if (adp<=tol){
                    d=2*dm*tol/(1+tol); alpha=-2*(sqrt(2/(1+tol))-1);
                }
                else{
                    d=2*trans/(dm+dp);
                    if (adm<=adp)
                        alpha=2*sqrt(d/dm)-2;
                    else
                        alpha=2-2*sqrt(d/dp);
                }
            }
        }

        //nu function
        if (fabs(alpha)<=tol){ nua=0.5; numa=0.5; }
        else{
            nua=(log((2-alpha)/(2+alpha))+2*alpha/(2-alpha))/(alpha*alpha);
            numa=(log((2+alpha)/(2-alpha))-2*alpha/(2+alpha))/(alpha*alpha);
        }

        //extrapolate
        result[k][0]=u[j][k]+d*h*nua;
        result[k][1]=u[j][k]-d*h*numa;
    }
}


void eno(double **result, double  **u, int j, int dim, double h){
    double dm, dp, d, Dm, Dp, D, der, Der;
    int k;
    for(k=0; k<dim; k++){
        dp=u[j+1][k]-u[j][k];
        dm=u[j][k]-u[j-1][k];
        d=0.5*(dp+dm);
        Dm=u[j][k]-2*u[j-1][k]+u[j-2][k];
        Dp=u[j+2][k]-2*u[j+1][k]+u[j][k];
        D=u[j+1][k]-2*u[j][k]+u[j-1][k];

        int choose;    

        //selection
        if (fabs(dm)<=fabs(dp)){
            if (fabs(Dm)<=fabs(D))
                choose=-1;
            else
                choose=0;
        }
        else{
            if (fabs(D)<=fabs(Dp))
                choose=0;
            else 
                choose=1;
        }

        //parabola parameters
        if (choose==-1){der=dm; Der=Dm;}
        if (choose==0){der=d; Der=D;}
        if (choose==1){der=dp; Der=Dp;}

        //extrapolate
        result[k][0]=u[j][k]+der*0.5+Der*(-1.0/24-choose*0.25+0.125);
        result[k][1]=u[j][k]-der*0.5+Der*(-1.0/24+choose*0.25+0.125);
    }
}


