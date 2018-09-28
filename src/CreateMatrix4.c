#include "inc.h"

float  CreateMatrix4(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,weight,data,cte)
float  **matrix,*g,*ge,*p,*pe;
int      setpro,setmat,checkpro;
char     seq[MAXP][MAXQ][MAXS];
profile *template,*target;
float    weight[MAXM],data[2][MAXM][MAXE][MAXE],cte[MAXP];
{

   int        i,j,k,n,ii,jj,lena,lenb;
   int        jtgt,itmp;
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
      for (ii=0;ii<template[n].size;ii++){
       for (jj=0;jj<target[n].size;jj++){
           jtgt=target[n].sequence[jj].element[j];
	   itmp=template[n].sequence[ii].element[i];
	   if ( (jtgt > 1 && itmp >1) || (jtgt <=1 && itmp <= 1) ){
            matrix[j][i] +=  cte[n-1]* (1.0/template[n].size) * (1.0/target[n].size) * data[1][n-1][jtgt][itmp];
	   }else{
            matrix[j][i] -=  cte[n-1]* (1.0/template[n].size) * (1.0/target[n].size) * data[1][n-1][jtgt][itmp];
           }
       }
      }
     }
    }
    for (k=0;k<setmat;k++){
     for (ii=0;ii<template[0].size;ii++){
      for (jj=0;jj<target[0].size;jj++){
          jtgt=target[0].sequence[jj].element[j];
          itmp=template[0].sequence[ii].element[i];
	   if ( (jtgt > 1 && itmp >1) || (jtgt <=1 && itmp <= 1) ){
            matrix[j][i] += ( 1.0 - ctp ) * weight[k] * (1.0/template[0].size) * (1.0/target[0].size) * data[0][k][jtgt][itmp];
	   }else{
            matrix[j][i] -= ( 1.0 - ctp ) * weight[k] * (1.0/template[0].size) * (1.0/target[0].size) * data[0][k][jtgt][itmp];
           }
      }
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
     for (ii=0;ii<template[n].size;ii++){
      itmp=template[n].sequence[ii].element[i];
      g[i]      +=  cte[n-1]*(1.0/template[n].size) * data[1][n-1][0][itmp];
      ge[i]     +=  cte[n-1]*(1.0/template[n].size) * data[1][n-1][1][itmp];
     }
     }
    }
    for (k=0;k<setmat;k++){
     for (ii=0;ii<template[0].size;ii++){
      itmp=template[0].sequence[ii].element[i];
      g[i]      += ( 1.0 - ctp ) * weight[k] * (1.0/template[0].size) * data[0][k][0][itmp];
      ge[i]     += ( 1.0 - ctp ) * weight[k] * (1.0/template[0].size) * data[0][k][1][itmp];
     }
    }
   }
   for (j=0;j<lenb;j++){
    p[j] =0.0;
    pe[j]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
     for (jj=0;jj<target[n].size;jj++){
      jtgt=target[n].sequence[jj].element[j];
      p[j]      +=  cte[n-1]*(1.0/target[n].size) *data[1][n-1][0][jtgt];
      pe[j]     +=  cte[n-1]*(1.0/target[n].size) *data[1][n-1][1][jtgt];
     }
     }
    }
    for (k=0;k<setmat;k++){
    for (jj=0;jj<target[0].size;jj++){
      jtgt=target[0].sequence[jj].element[j];
      p[j]      += ( 1.0 - ctp ) * weight[k] * (1.0/target[0].size) *data[0][k][0][jtgt];
      pe[j]     += ( 1.0 - ctp ) * weight[k] * (1.0/target[0].size) *data[0][k][1][jtgt];
    }
    }
   }
   return  media;
 }
