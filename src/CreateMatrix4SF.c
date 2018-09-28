#include "inc.h"

float  CreateMatrix4SF(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,weight,data,cte,lambda1,tolerance)
float  **matrix,*g,*ge,*p,*pe,*lambda1,*tolerance;
int      setpro,setmat,checkpro;
char     seq[MAXP][MAXQ][MAXS];
profile *template,*target;
psic   **PSICtemplate,**PSICtarget;
float    weight[MAXM],data[2][MAXM][MAXE][MAXE],cte[MAXP];
{

   int        i,j,k,n,ii,jj,lena,lenb;
   float      ctp,maximum,minimum,media,factor;
   

   lena=template[0].main.length;
   lenb=target[0].main.length;
   ctp=0.0;
   minimum=MAXL;
   maximum=0.0;


   factor=0.0;
   for (k=0;k<setpro-1;k++){ ctp += cte[k]; }
   if (ctp>1.0) ctp=1.0;


   for (i=0;i<lena;i++){
   for (j=0;j<lenb;j++){
    matrix[j][i]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
      for (ii=2;ii<MAXE;ii++){
       for (jj=2;jj<MAXE;jj++){
	   if ( (jj > 1 && ii >1) || (jj <=1 && ii <= 1) ){
            matrix[j][i] +=  cte[n-1]* PSICtemplate[n][i].frequency[ii] * PSICtarget[n][j].frequency[jj] * data[1][n-1][jj][ii]*lambda1[n];
	   }else{
            matrix[j][i] -=  cte[n-1]* PSICtemplate[n][i].frequency[ii] * PSICtarget[n][j].frequency[jj] * data[1][n-1][jj][ii]*lambda1[n];
           }
       }
      }
      factor+=cte[n-1]*lambda1[n];
     }
    }
    for (k=0;k<setmat;k++){
     n=k+setpro;
     if (ctp==1)break;
     for (ii=2;ii<MAXE;ii++){
      for (jj=2;jj<MAXE;jj++){
	   if ( (jj > 1 && ii >1) || (jj <=1 && ii <= 1) ){
            matrix[j][i] += ( 1.0 - ctp ) * weight[k] *  PSICtemplate[0][i].frequency[ii] * PSICtarget[0][j].frequency[jj] * data[0][k][jj][ii]*lambda1[n];
	   }else{
            matrix[j][i] -= ( 1.0 - ctp ) * weight[k] *  PSICtemplate[0][i].frequency[ii] * PSICtarget[0][j].frequency[jj] * data[0][k][jj][ii]*lambda1[n];
           }
      }
     }
     factor+=( 1.0 - ctp ) * weight[k] *lambda1[n];
    }
    if (matrix[j][i]>0 && minimum>matrix[j][i]){minimum=matrix[j][i];}
    if (matrix[j][i]>maximum){maximum=matrix[j][i];}
    }}



/*
   printf("             ");
   for (j=0;j<lenb;j++){printf("%8d:",j);}
   printf("\n");
   for(i=0;i<lena;i++){
   printf("ELEMENT%5d:",i);
   for (j=0;j<lenb;j++){
    if (matrix[j][i]>0 && minimum>matrix[j][i]){minimum=matrix[j][i];}
    if (matrix[j][i]>maximum){maximum=matrix[j][i];}
    printf("%8.1e;",matrix[j][i]);
   }printf("\n");}
*/

   if (minimum==MAXL)minimum=0.0;

   media=(maximum+minimum)/2.0;


   for (i=0;i<lena;i++){
    g[i] =0.0;
    ge[i]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
     for (ii=0;ii<MAXE;ii++){
      g[i]      +=  cte[n-1]*PSICtemplate[n][i].frequency[ii] * data[1][n-1][0][ii]* lambda1[n];
      ge[i]     +=  cte[n-1]*PSICtemplate[n][i].frequency[ii] * data[1][n-1][1][ii]* lambda1[n];
     }
     }
    }
    for (k=0;k<setmat;k++){
     for (ii=0;ii<MAXE;ii++){
      g[i]      += ( 1.0 - ctp ) * weight[k] * PSICtemplate[0][i].frequency[ii] * data[0][k][0][ii]* lambda1[k+setpro];
      ge[i]     += ( 1.0 - ctp ) * weight[k] * PSICtemplate[0][i].frequency[ii] * data[0][k][1][ii]* lambda1[k+setpro];
     }
    }
   }
   for (j=0;j<lenb;j++){
    p[j] =0.0;
    pe[j]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
     for (jj=0;jj<MAXE;jj++){
      p[j]      +=  cte[n-1]*PSICtarget[n][j].frequency[jj] *data[1][n-1][0][jj]* lambda1[n];
      pe[j]     +=  cte[n-1]*PSICtarget[n][j].frequency[jj] *data[1][n-1][1][jj]* lambda1[n];
     }
     }
    }
    for (k=0;k<setmat;k++){
    for (jj=0;jj<MAXE;jj++){
      p[j]      += ( 1.0 - ctp ) * weight[k] * PSICtarget[0][j].frequency[jj] *data[0][k][0][jj]* lambda1[k+setpro];
      pe[j]     += ( 1.0 - ctp ) * weight[k] * PSICtarget[0][j].frequency[jj] *data[0][k][1][jj]* lambda1[k+setpro];
    }
    }
   }

   return  media;
 }
