#include "inc.h"

float  CreateMatrix4PICASSO(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance)
float  **matrix,*g,*ge,*p,*pe,*lambda1,*tolerance;
int      setpro,setmat,checkpro;
char     seq[MAXP][MAXQ][MAXS];
profile *template,*target;
psic   **PSICtemplate,**PSICtarget;
float  **qtemplate,**qtarget;
float    weight[MAXM],data[2][MAXM][MAXE][MAXE],cte[MAXP];
{

   int        i,j,k,n,ii,jj,lena,lenb,dimset;
   float      ctp,maximum,minimum,media,factor;
   float      *lambda2,*pb,*qb;
   float    **score2;

   float      s;
   float      lambda();
   float     *fvector();
   float    **Fmatrix();
   void       free_Fmatrix();
   void       free_fvector();



   lena=template[0].main.length;
   lenb=target[0].main.length;


 
   ctp=0.0;
   minimum=MAXL;
   maximum=0.0;
   dimset=setpro+setmat;

   score2=Fmatrix(0,lena,0,lenb);
   lambda2=fvector(0,dimset);
   pb=fvector(0,lena);
   qb=fvector(0,lenb);
   

   for (i=0;i<lena;i++){for (j=0;j<lenb;j++){matrix[j][i]=0.0;}}

   for (i=0;i<setpro+setmat;i++){lambda2[i]=0.0;}
   

   lambda2[0]=1.0; 
   for (n=1;n<setpro;n++){
   if (cte[n-1]!=0){
    for (i=0;i<lena;i++){
      pb[i]=1.0/lena;
      for (j=0;j<lenb;j++){
        qb[j]=1.0/lenb;
        s=0.0;
        for (ii=2;ii<MAXE;ii++){
          if (qtemplate[n][ii]>0.0 && qtarget[n][ii] > 0.0){
             if (PSICtemplate[n][i].targetQ[ii]>0.0) s+=PSICtarget[n][j].frequency[ii]*log(PSICtemplate[n][i].targetQ[ii]/qtemplate[n][ii]);
             if (PSICtarget[n][j].targetQ[ii] >0.0)  s+=PSICtemplate[n][i].frequency[ii]*log(PSICtarget[n][j].targetQ[ii]/qtarget[n][ii]);
          }
        }
        score2[i][j] =  s ;
      }
    }
    if (lena>0 && lenb>0) lambda2[n]=lambda(pb,qb,score2,lena,lenb);
    printf("Lambda PICASSO[%d] = %e \n",n,lambda2[n]);
    if (lambda2[n]==0.0)cte[n-1]=0.0;
   }}

   for (n=setpro;n<setpro+setmat;n++){
   k=n-setpro;
   if (weight[k]!=0){
    for (i=0;i<lena;i++){
     pb[i]=1.0/lena;
     for (j=0;j<lenb;j++){
        qb[j]=1.0/lenb;
        s=0.0;
        for (ii=2;ii<MAXE;ii++){
          if (qtemplate[0][ii]>0.0 && qtarget[0][ii] > 0.0){
             if (PSICtemplate[n][i].targetQ[ii]>0.0) s+=PSICtarget[0][j].frequency[ii]*log(PSICtemplate[n][i].targetQ[ii]/qtemplate[0][ii]);
             if (PSICtarget[n][j].targetQ[ii] >0.0)  s+=PSICtemplate[0][i].frequency[ii]*log(PSICtarget[n][j].targetQ[ii]/qtarget[0][ii]);
          }
        }
        score2[i][j] =  s ;
      }
    }
    if (lena>0 && lenb>0) lambda2[n]=lambda(pb,qb,score2,lena,lenb);
    printf("Lambda PICASSO[%d] = %e \n",n,lambda2[n]);
    if (lambda2[n]==0.0)weight[k]=0.0;
   }}

   factor=0.0;
   for (k=0;k<setpro-1;k++){ ctp += cte[k]; }
   if (ctp>1.0) ctp=1.0;

   for (n=1;n<setpro;n++){
    for (i=0;i<lena;i++){
      for (j=0;j<lenb;j++){
        s=0.0;
        for (ii=2;ii<MAXE;ii++){
          if (qtemplate[n][ii]>0.0 && qtarget[n][ii] > 0.0){
             if (PSICtemplate[n][i].targetQ[ii]>0.0) s+=PSICtarget[n][j].frequency[ii]*log(PSICtemplate[n][i].targetQ[ii]/qtemplate[n][ii]);
             if (PSICtarget[n][j].targetQ[ii] >0.0)  s+=PSICtemplate[n][i].frequency[ii]*log(PSICtarget[n][j].targetQ[ii]/qtarget[n][ii]);
          }
        }
        matrix[j][i] +=  cte[n-1] * s * lambda2[n]  ;
      }
    }
    factor += cte[n-1] *lambda2[n];
    if (cte[n-1] > 0) printf("Multiplicative Factor[%d]= %e\n",n,cte[n-1] *lambda2[n]);
   }
   for (n=setpro;n<setpro+setmat;n++){
    k=n-setpro;
    if (ctp==1) break;
    for (i=0;i<lena;i++){
     for (j=0;j<lenb;j++){
        s=0.0;
        for (ii=2;ii<MAXE;ii++){
          if (qtemplate[0][ii]>0.0 && qtarget[0][ii] > 0.0){
             if (PSICtemplate[n][i].targetQ[ii]>0.0) s+=PSICtarget[0][j].frequency[ii]*log(PSICtemplate[n][i].targetQ[ii]/qtemplate[0][ii]);
             if (PSICtarget[n][j].targetQ[ii] >0.0)  s+=PSICtemplate[0][i].frequency[ii]*log(PSICtarget[n][j].targetQ[ii]/qtarget[0][ii]);
          }
        }
        matrix[j][i] +=  ( 1.0 - ctp ) * weight[k] * s * lambda2[n] ;
      }
    }
    factor +=  ( 1.0 - ctp ) * weight[k] *lambda2[n];
    if (( 1.0 - ctp ) * weight[k] > 0) printf("Multiplicative Factor[%d]= %e\n",n,( 1.0 - ctp ) * weight[k] *lambda2[n]);
   }

/*
   printf("             ");
   for (j=0;j<lenb;j++){printf("%8d-",j);}
   printf("\n");

   for(i=0;i<lena;i++){
   printf("ELEMENT %5d",i);
   for (j=0;j<lenb;j++){
    if (matrix[j][i]>0 && minimum>matrix[j][i]){minimum=matrix[j][i];}
    if (matrix[j][i]>maximum){maximum=matrix[j][i];}
    printf("%8.1e;",matrix[j][i]);
   }printf("\n");}
*/


   for(i=0;i<lena;i++){
   for (j=0;j<lenb;j++){
    if (matrix[j][i]>0 && minimum>matrix[j][i]){minimum=matrix[j][i];}
    if (matrix[j][i]>maximum){maximum=matrix[j][i];}
   }}


   if (minimum==MAXL)minimum=0.0;
   media=(maximum+minimum)/2.0;


   for (i=0;i<lena;i++){
    g[i] =0.0;
    ge[i]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
     for (ii=2;ii<MAXE;ii++){
      g[i]      +=  cte[n-1]*PSICtemplate[n][i].frequency[ii] * data[1][n-1][0][ii]* lambda1[n];
      ge[i]     +=  cte[n-1]*PSICtemplate[n][i].frequency[ii] * data[1][n-1][1][ii]* lambda1[n];
     }
     }
    }
    
    for (k=0;k<setmat;k++){
     for (ii=2;ii<MAXE;ii++){
      g[i]      += ( 1.0 - ctp ) * weight[k] * PSICtemplate[0][i].frequency[ii] * data[0][k][0][ii]* lambda1[k+setpro];
      ge[i]     += ( 1.0 - ctp ) * weight[k] * PSICtemplate[0][i].frequency[ii] * data[0][k][1][ii]* lambda1[k+setpro];
     }
    }
    if (g[i]  > fabs(maximum) ) g[i] =fabs(maximum);
    if (ge[i] > fabs(media/2) ) ge[i]=fabs(media/2);
   }
   for (j=0;j<lenb;j++){
    p[j] =0.0;
    pe[j]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
     for (jj=2;jj<MAXE;jj++){
      p[j]      +=  cte[n-1]*PSICtarget[n][j].frequency[jj] *data[1][n-1][0][jj]* lambda1[n] ;
      pe[j]     +=  cte[n-1]*PSICtarget[n][j].frequency[jj] *data[1][n-1][1][jj]* lambda1[n] ;
     }
     }
    }
    for (k=0;k<setmat;k++){
    for (jj=2;jj<MAXE;jj++){
      p[j]      += ( 1.0 - ctp ) * weight[k] * PSICtarget[0][j].frequency[jj] *data[0][k][0][jj]* lambda1[k+setpro];
      pe[j]     += ( 1.0 - ctp ) * weight[k] * PSICtarget[0][j].frequency[jj] *data[0][k][1][jj]* lambda1[k+setpro];
    }
    }
    if (p[i]  > fabs(maximum) ) p[i] =fabs(maximum);
    if (pe[i] > fabs(media/2) ) pe[i]=fabs(media/2);

   }

   free_fvector(lambda2,0,dimset);
   free_fvector(pb,0,lena);
   free_fvector(qb,0,lenb);
   free_Fmatrix(score2,0,lena,0,lenb);


   return  media;
 }


