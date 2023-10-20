#include "highOrder.h"

void reconstruct (int method, double **u, int j, int dim, double h, double *right, double *left){
    int i;
    
    //declare matrix to store output of extrapolation
    double **result;
    if ((result=(double**) malloc(dim*sizeof(double*)))==NULL) {printf("Error allocating memory for reconstruction\n"); exit(42);}
    for (i=0; i<dim; i++) if((result[i]=(double*) malloc(2*sizeof(double)))==NULL) {printf("Error allocating memory\n"); exit(42);}
    //reconstruct and extrapolate to interfaces left and right
    if (method==1){
        phm(result, u, j, dim, h);
    }
    if (method==2){
        eno(result, u, j, dim, h);
    }

    for (i=0;i<dim;i++){
        right[i]=result[i][0];
        left[i]=result[i][1];
    }
    //free used memory
    for (i=0; i<dim; i++){
        free(result[i]);
    }
    free(result);
}