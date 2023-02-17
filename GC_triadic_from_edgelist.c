#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define N 10000
#define cm 13
#define Nrunmax 100
#define gamma 2.5
#define m    1
#define kc   3
//#define p0   0.6
#define cn   1.7
#define cp   10

#define err_tol 0.01
int *vis1,*size_cluster1,**knn1,*k1,c1,c2,c3,*occ;
int **y;
float p0;


/**************************************************************************
 Recurrence is a subrutine for calulating the giant component in link percolation
 **************************************************************************/
int Recurrence( int i , int cluster_size, int ncluster){
	int j, n3,aus_cluster_size;
	
	cluster_size++;		
		vis1[i]=ncluster;
		for(n3=0;n3<k1[i];n3++){
			j=knn1[i][n3];
			if((y[i][j]==1)&&(vis1[j]==0)){
				aus_cluster_size=Recurrence(j, cluster_size, ncluster);
				cluster_size=aus_cluster_size;
			}
	
		}
	
		return cluster_size;
}

int main(int argc, char** argv){
	int i,j,n,it,ncluster1,ncluster2,ncluster3, GMCC, cluster_size,m1,m2,m3,m1_aus,m2_aus,m3_aus,c1_aus,c2_aus,c3_aus,*sigma,nrun,nc,np,nc2,niter,j2,ni2;
	int s1,s2,Nc1,Nc2,Nc3,aus,aus3,**adj1,**adj2,**adj3,N0,**s,nrun2,**denR,*T,**knr,**knp,***knn_np,***knn_nr;
	float p,f,**x,*Sm,**n1,**n2,MGC1,MGC2,MGC3,MGC4,nsum1,nsumold1,nsum2,nsumold2,aus1,aus2,sigma11,sigma10,*sigma11m,*sigma10m,*kx,Norm1,**GCA,rho,r,*Tx;
    int ausi,GC;
    float pc,kav,kappa,kappaT,**pL;
    double x2;
	   char filename[60],string1[50],string2[50];
	
	FILE *gp2,*gp,*fp,*gp3;
  

  /**************************************************************************
   open file for output 
   finemane GC
   at the end of teh program the file will contain 
   two columns: p GC GC2
   
   **************************************************************************/

		srand48(time(NULL));
	//filename="british_edgelist.txt";
	
		
	
	N0=N;
	//printf("%d\n",N);
	
	Tx=(float*)calloc(N,sizeof(float));
	vis1=(int*)calloc(N,sizeof(int));
    kx=(float*)calloc(N,sizeof(float));

	occ=(int*)calloc(N,sizeof(int));
	k1=(int*)calloc(N,sizeof(int));
	x=(float**)calloc(N,sizeof(float*));
    y=(int**)calloc(N,sizeof(int*));
	T=(int*)calloc(N,sizeof(int));
	
    knr=(int**)calloc(N,sizeof(int*));
    knp=(int**)calloc(N,sizeof(int*));
    knn_nr=(int***)calloc(N,sizeof(int**));
    knn_np=(int***)calloc(N,sizeof(int**));
	knn1=(int**)calloc(N,sizeof(int*));
	pL=(float**)calloc(N,sizeof(float*));

	adj1=(int**)calloc(N,sizeof(int*));

        GCA=(float**)calloc(100,sizeof(float*));
		for(i=0;i<N;i++){
			knn1[i]=(int*)calloc(N,sizeof(int));
            knr[i]=(int*)calloc(N,sizeof(int));
            knp[i]=(int*)calloc(N,sizeof(int));
            knn_nr[i]=(int**)calloc(N,sizeof(int*));
             knn_np[i]=(int**)calloc(N,sizeof(int*));
			adj1[i]=(int*)calloc(N,sizeof(int));
			x[i]=(float*)calloc(N,sizeof(float));
            y[i]=(int*)calloc(N,sizeof(int));
       

            pL[i]=(float*)calloc(N,sizeof(float));
            if(i<100)
            GCA[i]=(float*)calloc(100,sizeof(float));
            
		}
	size_cluster1=(int*)calloc(N,sizeof(int));



	
	for(i=0;i<N;i++){
		k1[i]=0;
	}

    printf("ciao 0\n");

    sprintf(filename,"edge_list.txt");
    fp=fopen(filename,"r");
    
    while(!feof(fp)){
        fscanf(fp,"%d %d ", &i ,&j);
				k1[i]++;
				k1[j]++;
				knn1[i][k1[i]-1]=j;
				knn1[j][k1[j]-1]=i;
                adj1[i][j]=1;
                adj1[j][i]=1;
            knn_nr[i][j]=(int*)calloc(N,sizeof(int));
            knn_np[i][j]=(int*)calloc(N,sizeof(int));
            knn_nr[j][i]=(int*)calloc(N,sizeof(int));
            knn_np[j][i]=(int*)calloc(N,sizeof(int));
            knr[i][j]=0;
            knp[i][j]=0;
            knr[j][i]=0;
            knp[j][i]=0;
    }
    printf("ciao 2\n");
    sprintf(filename,"triadic_edge_list_plus.txt");

    fp=fopen(filename,"r");
    printf("ciao 2b\n");
    while(!feof(fp)){
        fscanf(fp,"%d %d %d ", &i ,&j, &n);
     //   printf("%d %d %d \n", i,j,n);
        knp[i][j]++;
        knn_np[i][j][knp[i][j]-1]=n;
        knp[j][i]++;
        knn_np[j][i][knp[j][i]-1]=n;
        
    }
        
    printf("ciao 3\n");
    sprintf(filename,"triadic_edge_list_minus.txt");
    fp=fopen(filename,"r");
    
    while(!feof(fp)){
        fscanf(fp,"%d %d %d", &i ,&j, &n);
                knr[i][j]++;
                knn_nr[i][j][knr[i][j]-1]=n;
                knr[j][i]++;
                knn_nr[j][i][knr[j][i]-1]=n;
                
    }
        
		
     
     fclose(fp);
			
    printf("ciao\n");
    sprintf(filename,"GC.txt");
    gp2=fopen(filename,"w");
    for(p0=1;p0>0.01;p0=p0-0.01){
	
	for (niter=0;niter<Nrunmax;niter++){
     printf("niter=%d\n",niter);
        
        for(i=0;i<N;i++){
            
            for(n=0;n<k1[i];n++){
                j=knn1[i][n];
                
               // if(j<i){
                    x[i][j]=drand48();
                    x[j][i]=x[i][j];
                //}
            }
        }
    

        nc=0;
        
        for(i=0;i<N;i++){
            for(n=0;n<k1[i];n++){
                 j=knn1[i][n];
                if(niter>1){
                aus=1.;
                aus2=1.;
                   
                for(ni2=0;ni2<knr[i][j];ni2++){
                    aus=aus*(float)(1-occ[knn_nr[i][j][ni2]]);
                }
                for(ni2=0;ni2<knp[i][j];ni2++){
                aus2=aus2*(float)(1-occ[knn_np[i][j][ni2]]);
                }
                }
            if(niter==0){
                
                aus=1;aus2=0;}
               // aus=1;aus2=0;
        
            pL[i][j]=p0*aus*(1.-aus2);
            }
        }
    
        GC=0.;
        for(i=0;i<N;i++){
            vis1[i]=0;
        }

        for(i=0;i<N;i++){
            for(n=0;n<k1[i];n++){
                j=knn1[i][n];
                y[i][j]=0;
                y[j][i]=0;
                if((x[i][j]<pL[i][j])){
                y[i][j]=1;
                y[j][i]=1;
                } 
            }
        }
        
    
		ncluster1=0;		
		for(i=0;i<N;i++){
			vis1[i]=0;
		}
		m1=0;
		ncluster1=0;
            c1=0;
		for(n=0;n<N;n++){
			if(vis1[n]==0){
				cluster_size=0;
				ncluster1++;
				cluster_size=Recurrence(n, cluster_size, ncluster1);
				size_cluster1[ncluster1]=cluster_size;
				if(cluster_size>m1){m1=cluster_size;
					c1=ncluster1;}				
			}			
		}
        
		
		Nc1=c1;
			//m2=0;
            GC=0.;
			for(i=0;i<N;i++){
                occ[i]=0;
				if(vis1[i]==Nc1){
                    occ[i]=1;
					GC++;
				//	m2++;
				}
                
                
            }
//        denR[nc][GC]++;
 //       GCA[nc][nc]+=(float)GC;
        printf("%f,\n",(float)GC);
        if(niter>20){
            fprintf(gp2,"%f %d %f \n",p0,niter,((float)GC)/(float)N);
    }
    }
    }
    //  fclose(gp3);
    return 0;
}
    
		
