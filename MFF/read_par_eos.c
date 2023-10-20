#include "read_par_eos.h"

char *strip_copy(char const *s)
{
    char *buf = malloc(1 + strlen(s));
    if (buf)
    {
        char *p = buf;
        char const *q;
        int n;
        for (q = s; *q; q += n + strspn(q+n, "\n"))
        {
            n = strcspn(q, "\n");
            strncpy(p, q, n);
            p += n;
        }
        *p++ = '\0';
        buf = realloc(buf, p - buf);
    }
    return buf;
}

void read_par_eos(struct PP *eosPP, struct TASPP *eosTASPP, struct Tab *eosTab, char *parfile, int *usePP, int *useTab, int *useTASPP){
    FILE *fp = fopen(parfile, "r");
    if(fp==NULL){printf("Error opening eos.par\n"); exit(1);}
    const unsigned MAX_LENGTH = 256;
    char buffer[MAX_LENGTH], line2[MAX_LENGTH];
    int count=0;
    char *token;
    int whichEoS; char *path; int n_PT; double *denter, *xmenter, *rhoPT, **paramsenter, *Kfitted, *Gfitted, extra_cut, rc;

    fgets(buffer, MAX_LENGTH, fp); //skip header
    
    //Read type of Eos
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    if (strstr(line2, "TASPP") != NULL) {
        (*useTASPP)=1;
    }
    else if (strstr(line2, "PP") != NULL) {
        (*usePP)=1;
    }
    else if(strstr(line2, "TAB") != NULL) {
        (*useTab)=1;
    }

    //Go to next field
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    path=strip_copy(line2);
    
    //If EoS is PP
    if (*usePP){
        whichEoS=stringToEoS(path);
        preparePP(eosPP, whichEoS); 
    }
    //If EoS is Tabulated
    if (*useTab){ 
        prepareTab(path, eosTab);
    }
    //If EoS is TASPP
    if (*useTASPP){ 
        n_PT=atoi(path);
        denter=(double*)malloc(n_PT*sizeof(double)); xmenter=(double*)malloc(n_PT*sizeof(double));
        rhoPT=(double*)malloc(2*n_PT*sizeof(double)); paramsenter=(double**)malloc(n_PT*sizeof(double*));
        for(int i=0; i<n_PT; i++) paramsenter[i]=(double*)malloc(4*sizeof(double));

        fgets(buffer, MAX_LENGTH, fp); //skip line
        fgets(buffer, MAX_LENGTH, fp); //skip header
        fgets(buffer, MAX_LENGTH, fp); //skip header
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        token = strtok(line2, "\t");
        count=0;
        while( token != NULL ) {
            denter[count]=strtod(token,NULL);
            token = strtok(NULL, "\t");
            count++;
        }
        count=0;
        fgets(buffer, MAX_LENGTH, fp); //skip line
        fgets(buffer, MAX_LENGTH, fp); //skip header
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        token = strtok(line2, "\t");
        while( token != NULL ) {
            xmenter[count]=strtod(token,NULL);
            token = strtok(NULL, "\t");
            count++;
        } 
        count=0;
        fgets(buffer, MAX_LENGTH, fp); //skip line
        fgets(buffer, MAX_LENGTH, fp); //skip header
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        token = strtok(line2, "\t");
        while( token != NULL ) {
            rhoPT[count]=strtod(token,NULL);
            token = strtok(NULL, "\t");
            count++;
        }
        count=0;
        fgets(buffer, MAX_LENGTH, fp); //skip line
        fgets(buffer, MAX_LENGTH, fp); //skip header
        for(int j=0; j<n_PT;j++){
            fgets(buffer, MAX_LENGTH, fp);
            memcpy(line2, buffer, MAX_LENGTH);
            token = strtok(line2, "\t");
            while( token != NULL ) {
                paramsenter[j][count]=strtod(token,NULL);
                token = strtok(NULL, "\t");
                count++;
            }
            count=0;
        }
        fgets(buffer, MAX_LENGTH, fp); //skip line
        fgets(buffer, MAX_LENGTH, fp); //skip header
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        token = strtok(line2, "\t");
        rc=strtod(token,NULL);
        token = strtok(NULL, "\t");
        extra_cut=strtod(token,NULL);

        Kfitted=(double*)malloc((n_PT+1+SIGN0(extra_cut))*sizeof(double));
        Gfitted=(double*)malloc((n_PT+1+SIGN0(extra_cut))*sizeof(double));

        fgets(buffer, MAX_LENGTH, fp); //skip line
        fgets(buffer, MAX_LENGTH, fp); //skip header
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        count=0;
        token = strtok(line2, "\t");
        while( token != NULL ) {
            Kfitted[count]=strtod(token,NULL);
            token = strtok(NULL, "\t");
            count++;
        }

        fgets(buffer, MAX_LENGTH, fp); //skip line
        fgets(buffer, MAX_LENGTH, fp); //skip header
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        count=0;
        token = strtok(line2, "\t");
        while( token != NULL ) {
            Gfitted[count]=strtod(token,NULL);
            token = strtok(NULL, "\t");
            count++;
        }
        prepareTASPP(rhoPT, paramsenter, denter, xmenter, n_PT, eosTASPP, extra_cut, rc, Kfitted, Gfitted);
    }
    fclose(fp);
}




    

    

    
