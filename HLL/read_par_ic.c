#include "read_par_ic.h"

void read_par_ic(char *parfile, char **filename, double *dInfo, int *iInfo){
    FILE *fp = fopen(parfile, "r");
    if(fp==NULL){printf("Error opening ic.par\n"); exit(1);}
    const unsigned MAX_LENGTH = 256;
    char buffer[MAX_LENGTH], line2[MAX_LENGTH];

    fgets(buffer, MAX_LENGTH, fp); //skip header
    
    //Initial conditions
    for(int k=0; k<4; k++){
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        dInfo[k]=strtod(line2,NULL);
    }
    

    //spatial domain
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    for(int k=4; k<6; k++){
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        dInfo[k]=strtod(line2,NULL);
    }

    //spatial discretization
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header    
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    iInfo[0]=atoi(line2);

    //cfl condition
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header    
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    dInfo[6]=strtod(line2,NULL);

    //final time
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header    
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    dInfo[7]=strtod(line2,NULL);

    //output file name
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH); 
    (*filename)=strip_copy(line2); 

    //allowUpwind and interpolation method
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header    
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    iInfo[1]=atoi(line2);
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    iInfo[2]=atoi(line2);
    
}


