/**************************************************************************************************
 * If you use this code, please cite:
 *
 * 1. 
 * G. Bianconi and O.T. Courtney
 * "Generalized network structures: the configuration model and the canonical ensemble of
 * simplicial complexes"
 * Phys. Rev. E 93, 062311 (2016)
 * 
 * 2.  
 * Hanlin Sun, Filippo Radicchi, Jürgen Kurths and Ginestra Bianconi
 * "The dynamic nature of percolation on networks with triadic interactions"
 * arXiv:2204.13067

***************************************************************************************************
 * Code that structural networks with scale-free degree distribution and regulatory networks
 * with Poisson degree distribution.
 *
 * The option to use a Poisson strucutural network has also been included.
 * The necessary code may be found in comments at the relevant points.
 *
 * This code uses:
 * N  Number of nodes in the structural network
 * m  The minimum of the scale-free distribution
 * gamma2  Exponent of the scale-free distribution
 * lambda  Expected value of thes structural Poisson distribution (commented-out)
 * Avoid  Whether or not 'back-tracking' is allowed when illegal matchings are proposed
 * (Avoid==1 allowed, Avoid==0 not allowed)
 * NX  Maximum number of 'back-tracks' before matching process restarts from an unmatched network
 * K maximum allowed degree
 * cp Expected value of positive regulatory Poisson distribution
 * cn Expected value of negative regulatory Poisson distribution
 *************************************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define N 10000
#define m 4
#define gamma2 2.5
/* #define lambda 10 */
#define Avoid 1
#define NX 90
#define K 100
#define cp 10.
#define cn 2.8

int *kgi,*kg,***tri;

/*************************************************************************************************/
/* Randomly select an unmatched stub. Choose takes as its input a random number between 0 and the
total number of stubs and gives as its output the index of the node of the selected stub */
int Choose(double x){
    int i1,i;
    for (i=0;i<N;i++){
        x-=kgi[i];
        if (x<0){
            i1=i;
            break;
        }
    }
    return(i1);
}
/*************************************************************************************************/

int main(int argc, char** argv){
    int i,j,nrun,j2,i1,i2,i3,naus,*knng,*pkg,*k,**l,*pk,*knn,n,**a,*Ck;
    double xaus, x;
    char filec[60];

    FILE *fp,*gp,*gp2,*gp3;

    srand48(time(NULL));
    kgi=(int*)calloc(N,sizeof(int));
    kg=(int*)calloc(N,sizeof(int));
    k=(int*)calloc(N,sizeof(int));
    a=(int**)calloc(N,sizeof(int*));
    knng=(int*)calloc(N,sizeof(int));
    pkg=(int*)calloc(N,sizeof(int));
    knn=(int*)calloc(N,sizeof(int));
    pk=(int*)calloc(N,sizeof(int));
    Ck=(int*)calloc(N,sizeof(int));

    for(i=0;i<N;i++){
        a[i]=(int*)calloc(N,sizeof(int));
    }

    xaus=4;

    while(xaus>2){
    /***********************************************************************************************/
    /* Initialization */
        for(i=0;i<N;i++){
        /* Nodes are assigned desired generalized degree according to a scale-free distribution */
            kgi[i]=(int)(m*pow(drand48(),-1./(gamma2-1.)));
            /* kgi[i]= poisson(lambda); */
            while(kgi[i]>(K)){
            /* Desired generalized degrees are re-drawn if they exceed the maximum possible generalized degree of a node (natural cut-off) */
                kgi[i]=(int)(m*pow(drand48(),-1./(gamma2-1.)));
                /* kgi[i]= poisson(lambda); */
            }
            kg[i]=0;  /* Generalized degree of node i intially set to 0 */
            k[i]=0;  /* Degree of node i intially set to 0 */
            for(j=0;j<N;j++){
                a[i][j]=0;
            }
        }
        xaus=0;
        for(i=0;i<N;i++){
            xaus+=kgi[i];
        }
        naus=0; /* Back-track counter initially set to zero */
    /***********************************************************************************************/
    /* Stubs matched */
        while((xaus>3)&&(naus<1+Avoid*NX)){
            /* Randomly select two nodes proportional to the number of unmatched stubs they have remaining. */
            x=xaus*drand48();
            i1=Choose(x);
            kg[i1]++;
            kgi[i1]--;
            xaus--;

            x=xaus*drand48();
            i2=Choose(x);
            kg[i2]++;
            kgi[i2]--;
            xaus--;

            /* Check proposed matching is legal */
            if((i1!=i2)&&(a[i1][i2]==0)){
                /* Proposed matching legal. Create link */
                a[i1][i2]=1;
                a[i2][i1]=1;
            }
            else{
                /* Proposed matching illegal. Back-track and increment back-track counter by one */
                naus++;
                if(Avoid==1){
                    kg[i1]--;
                    kgi[i1]++;
                    kg[i2]--;
                    kgi[i2]++;
                }
            }
        }
    }
/*************************************************************************************************/
/* Degrees calculated */
    for (i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            if(a[i][j]>0){
                k[i]++;
                k[j]++;
            }
        }
    }
/*************************************************************************************************/
/* Print list of edges to file */
    
    gp=fopen("edge_list.txt","w");
    gp2=fopen("triadic_edge_list_plus.txt","w");
    gp3=fopen("triadic_edge_list_minus.txt","w");

/* Generate the positive and negative triadic regulation. Edge (i,j) is positively/negatively regulated
   by node n. */
    for (i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            if(a[i][j]==1){
                fprintf(gp,"%d %d\n",i,j);
                for(n=0;n<N;n++){
                    xaus=drand48();
                        if(xaus<(cp)/((float)N)){
                            fprintf(gp2,"%d %d %d\n",i,j,n);
                        }
                        if((xaus>(cp)/((float)N))&&(xaus<(cp+cn)/((float)N))){
                            fprintf(gp3,"%d %d %d\n",i,j,n);
                        }
                }
            }
        }
    }
        


/*************************************************************************************************/
    fclose(gp);
    fclose(gp2);
    fclose(gp3);

    return 0;
}


   


