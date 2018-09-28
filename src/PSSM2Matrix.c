#include  "inc.h"

float PSSM2Matrix(pssm,matrix,g,ge,p,pe,setpro,setmat,checkpro,template,target,qtemplate,PSICtemplate,PSICtarget,qtarget,weight,data,cte,lambda1,tolerance)
float  **pssm,**matrix,*g,*ge,*p,*pe,*lambda1,*tolerance;
int      setpro,setmat,checkpro;
profile *template,*target;
psic   **PSICtemplate,**PSICtarget;
float  **qtemplate,**qtarget;
float    weight[MAXM],data[2][MAXM][MAXE][MAXE],cte[MAXP];
{
   int        lena,lenb,i,j,n,ii,jj,k;
   float     *pb,*qb,**score2;
   float      lambda2,ctp,maximum,minimum,media,maxgap,maxgapex,mingap,mingapex,maximum2;
   float      lambda();
   float     *fvector();
   float    **Fmatrix();
   void       free_Fmatrix();
   void       free_fvector();


   printf("\n**************************\n");
   printf("\nMETHOD:\tMENTA \n");

   lena=template[0].main.length;
   lenb=target[0].main.length;
   ctp=0.0;
   minimum=mingap=mingapex=MAXL;
   maximum=maxgap=maxgapex=0.0;
   media=0.0;

   for (i=0;i<lena;i++){for (j=0;j<lenb;j++){matrix[j][i]=0.0;}}
   for (i=0;i<lena;i++){g[i] =0.0; ge[i]=0.0;}
   for (j=0;j<lenb;j++){p[j] =0.0; pe[j]=0.0;}

   score2=Fmatrix(0,lena,0,lenb);
   pb=fvector(0,lena);
   qb=fvector(0,lenb);
   
   for (k=0;k<setpro-1;k++){ ctp += cte[k]; }
   if (ctp>1.0) ctp=1.0; 

   lambda2=1.0; 
   for (i=0;i<lena;i++){
      pb[i]=1.0/lena;
      for (j=0;j<lenb;j++){
        qb[j]=1.0/lenb;
        score2[i][j] =  pssm[j+1][i+1] ;
        matrix[j][i] =  0.0;
      }
   }
   lambda2=lambda(pb,qb,score2,lena,lenb);

   for (n=1;n<setpro;n++){
   if (cte[n-1] > 0){
    for (i=1;i<=lena;i++){
    for (j=1;j<=lenb;j++){
     matrix[j-1][i-1] += cte[n-1] * pssm[j][i] * lambda2;
    }}
    if (cte[n-1] > 0) printf("Multiplicative Factor[%d]= %e\n",n,cte[n-1] *lambda2);
   }}
   for (n=setpro;n<setpro+setmat;n++){
   k=n-setpro;
   if (ctp==1) break;
   if (weight[k] > 0){
    for (i=1;i<=lena;i++){
    for (j=1;j<=lenb;j++){
     matrix[j-1][i-1] += ( 1.0 - ctp ) * weight[k] * pssm[j][i] * lambda2;
    }}
    if (( 1.0 - ctp ) * weight[k] > 0) printf("Multiplicative Factor[%d]= %e\n",n,( 1.0 - ctp ) * weight[k] *lambda2);
   }}


   tolerance[6]=lambda2*tolerance[6];

   for (i=0;i<lena;i++){
   for (j=0;j<lenb;j++){
    if (matrix[j][i]>0 && minimum>matrix[j][i]){minimum=matrix[j][i];}
    if (matrix[j][i]>maximum){maximum=matrix[j][i];}
    if (matrix[j][i]>maximum2 && matrix[j][i]<maximum){maximum2=matrix[j][i];}
    }}

   media=(maximum2+minimum)/2.0;

   printf("Lambda Score   =\t%e \n",lambda2);
   for (n=1;n<setpro;n++){printf("Lambda Gap(F)[%3d]=\t%e\n",n,cte[n-1]*lambda1[n]);}
   for (k=0;k<setmat;k++){printf("Lambda Gap(W)[%3d]=\t%e\n",k,( 1.0 - ctp ) * weight[k] * lambda1[k+setpro]);}
   printf("Maximum(1)     =\t%e\n",maximum);
   printf("Maximum(2)     =\t%e\n",maximum2);
   printf("Minimum        =\t%e\n",minimum);

   for (i=0;i<lena;i++){
    g[i] =0.0;
    ge[i]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
     for (ii=0;ii<MAXE;ii++){
      g[i]      +=  cte[n-1]*PSICtemplate[n][i].frequency[ii] * data[1][n-1][0][ii] * lambda1[n] ;
      ge[i]     +=  cte[n-1]*PSICtemplate[n][i].frequency[ii] * data[1][n-1][1][ii] * lambda1[n];
     }
     }
    }
    for (k=0;k<setmat;k++){
     for (ii=0;ii<MAXE;ii++){
      g[i]      += ( 1.0 - ctp ) * weight[k] * PSICtemplate[0][i].frequency[ii] * data[0][k][0][ii] * lambda1[k+setpro];
      ge[i]     += ( 1.0 - ctp ) * weight[k] * PSICtemplate[0][i].frequency[ii] * data[0][k][1][ii] * lambda1[k+setpro];
     }
    }
    if (g[i]  > fabs(maximum2) ) g[i] =fabs(maximum2);
    if (ge[i] > fabs(media/2) ) ge[i]=fabs(media/2);
    if (g[i]>maxgap)maxgap=g[i];
    if (ge[i]>maxgapex)maxgapex=ge[i];
    if (g[i]<mingap && g[i]>0.0 )mingap=g[i];
    if (ge[i]<mingapex && ge[i]>0.0 )mingapex=ge[i];
   }

   for (j=0;j<lenb;j++){
    p[j] =0.0;
    pe[j]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
     for (jj=0;jj<MAXE;jj++){
      p[j]      +=  cte[n-1]*PSICtarget[n][j].frequency[jj] *data[1][n-1][0][jj] * lambda1[n];
      pe[j]     +=  cte[n-1]*PSICtarget[n][j].frequency[jj] *data[1][n-1][1][jj] * lambda1[n];
     }
     }
    }
    for (k=0;k<setmat;k++){
    for (jj=0;jj<MAXE;jj++){
      p[j]      += ( 1.0 - ctp ) * weight[k] * PSICtarget[0][j].frequency[jj] *data[0][k][0][jj] * lambda1[k+setpro];
      pe[j]     += ( 1.0 - ctp ) * weight[k] * PSICtarget[0][j].frequency[jj] *data[0][k][1][jj] * lambda1[k+setpro];
    }
    }
    if (p[j]  > fabs(maximum2) ) p[j] =fabs(maximum2);
    if (pe[j] > fabs(media/2) ) pe[j]=fabs(media/2);
    if (p[i]>maxgap)maxgap=p[i];
    if (pe[i]>maxgapex)maxgapex=pe[i];
    if (p[i]<mingap && p[i]>0.0 )mingap=p[i];
    if (pe[i]<mingapex && pe[i]>0.0)mingapex=pe[i];
   }

   printf("Maximum Gap( & ext.) =\t%e\t%e\n",maxgap,maxgapex);
   printf("Minimum Gap( & ext.) =\t%e\t%e\n",mingap,mingapex);

   

   free_fvector(pb,0,lena);
   free_fvector(qb,0,lenb);
   free_Fmatrix(score2,0,lena,0,lenb);

/*   printf("Media   =\t%e\n",media); */

   return  media;
}
