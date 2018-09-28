#include "inc.h"

float  CreateMatrix4DPQ(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,weight,data,cte)
float  **matrix,*g,*ge,*p,*pe;
int      setpro,setmat,checkpro;
char     seq[MAXP][MAXQ][MAXS];
profile *template,*target;
psic   **PSICtemplate,**PSICtarget;
float    weight[MAXM],data[2][MAXM][MAXE][MAXE],cte[MAXP];
{

   int        i,j,k,n,ii,jj,lena,lenb;
   float      ctp,maximum,minimum,media;


   lena=template[0].main.length;
   lenb=target[0].main.length;
   ctp=0.0;
   minimum=MAXL;
   maximum=0.0;

   for (k=0;k<setpro-1;k++){ ctp += cte[k]; }

   for (i=0;i<lena;i++){
   for (j=0;j<lenb;j++){
    matrix[j][i]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
      for (ii=0;ii<MAXE;ii++){
            matrix[j][i] +=  cte[n-1]* PSICtemplate[n][i].targetQ[ii] * PSICtarget[n][j].targetQ[ii] ;
      }
     }
    }
    for (k=0;k<setmat;k++){
     for (ii=0;ii<MAXE;ii++){
            matrix[j][i] += ( 1.0 - ctp ) * weight[k] *  PSICtemplate[0][i].targetQ[ii] * PSICtarget[0][j].targetQ[ii] ;
     }
    }
    if (matrix[j][i]>0 && minimum>matrix[j][i]){minimum=matrix[j][i];}
    if (matrix[j][i]>maximum){maximum=matrix[j][i];}
    }}

   media=(maximum+minimum)/2.0;
   for (i=0;i<lena;i++){
    g[i] =0.0;
    ge[i]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
     for (ii=0;ii<MAXE;ii++){
      g[i]      +=  cte[n-1]*PSICtemplate[n][i].targetQ[ii] * data[1][n-1][0][ii];
      ge[i]     +=  cte[n-1]*PSICtemplate[n][i].targetQ[ii] * data[1][n-1][1][ii];
     }
     }
    }
    for (k=0;k<setmat;k++){
     for (ii=0;ii<MAXE;ii++){
      g[i]      += ( 1.0 - ctp ) * weight[k] * PSICtemplate[0][i].targetQ[ii] * data[0][k][0][ii];
      ge[i]     += ( 1.0 - ctp ) * weight[k] * PSICtemplate[0][i].targetQ[ii] * data[0][k][1][ii];
     }
    }
   }
   for (j=0;j<lenb;j++){
    p[j] =0.0;
    pe[j]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
     for (jj=0;jj<MAXE;jj++){
      p[j]      +=  cte[n-1]*PSICtarget[n][j].targetQ[jj] *data[1][n-1][0][jj];
      pe[j]     +=  cte[n-1]*PSICtarget[n][j].targetQ[jj] *data[1][n-1][1][jj];
     }
     }
    }
    for (k=0;k<setmat;k++){
    for (jj=0;jj<MAXE;jj++){
      p[j]      += ( 1.0 - ctp ) * weight[k] * PSICtarget[0][j].targetQ[jj] *data[0][k][0][jj];
      pe[j]     += ( 1.0 - ctp ) * weight[k] * PSICtarget[0][j].targetQ[jj] *data[0][k][1][jj];
    }
    }
   }
   return  media;
 }
