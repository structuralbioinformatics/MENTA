#include "inc.h"

float  CreateMatrix4DPF(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte)
float  **matrix,*g,*ge,*p,*pe;
int      setpro,setmat,checkpro;
char     seq[MAXP][MAXQ][MAXS];
profile *template,*target;
psic   **PSICtemplate,**PSICtarget;
float  **qtemplate,**qtarget;
float    weight[MAXM],data[2][MAXM][MAXE][MAXE],cte[MAXP];
{

   int        i,j,k,n,ii,jj,lena,lenb;
   float      ctp,maximum,minimum,media;
   float      *lambda1,*lambda2,*pb,*qb,s;
   float    **score1,**score2;
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
   score1=Fmatrix(0,MAXE,0,MAXE);
   score2=Fmatrix(0,lena,0,lenb);
   lambda1=fvector(0,setpro+setmat);
   lambda2=fvector(0,setpro+setmat);
   pb=fvector(0,MAXE+lena);
   qb=fvector(0,MAXE+lenb);

   for (k=0;k<setpro-1;k++){ ctp += cte[k]; }


   lambda1[0]=1.0;
   for (n=1;n<setpro;n++){
    pb[0]=qb[0]=pb[1]=qb[1]=0.0;
    for (ii=2;ii<MAXE;ii++){
       score1[0][ii]=score1[ii][0]=score1[1][ii]=score1[ii][1]=0.0;
       pb[ii]=qtemplate[n][ii];
       qb[ii]=qtarget[n][ii];
       for (jj=2;jj<MAXE;jj++){score1[ii][jj]=data[1][n-1][ii][jj];}
    }
    lambda1[n]=lambda(pb,qb,score1,MAXE,MAXE);
    printf("lambda %d = %f\n",n,lambda1[n]);
   }
   for (n=setpro;n<setpro+setmat;n++){
    k=n-setpro;
    pb[0]=qb[0]=pb[1]=qb[1]=0.0;
    for (ii=2;ii<MAXE;ii++){
       score1[0][ii]=score1[ii][0]=score1[1][ii]=score1[ii][1]=0.0;
       pb[ii]=qtemplate[0][ii];
       qb[ii]=qtarget[0][ii];
       for (jj=2;jj<MAXE;jj++){score1[ii][jj]=data[0][k][ii][jj];}
    }
    lambda1[n]=lambda(pb,qb,score1,MAXE,MAXE);
    printf("lambda %d = %f\n",n,lambda1[n]);
   }   
  
   lambda2[0]=1.0; 
   for (n=1;n<setpro;n++){
    pb[0]=qb[0]=pb[1]=qb[1]=0.0;
    for (i=0;i<lena;i++){
      pb[i]=1.0/lena;
      for (j=0;j<lenb;j++){
        qb[j]=1.0/lenb;
        s=0.0;
        for (ii=2;ii<MAXE;ii++){
             s+= (PSICtemplate[n][i].frequency[ii] * PSICtarget[n][j].frequency[ii]);
        }
        score2[i][j] =  s ;
      }
    }
    lambda2[n]=lambda(pb,qb,score2,lena,lenb);
    printf("lambda2 %d = %f\n",n,lambda2[n]);
   }
   for (n=setpro;n<setpro+setmat;n++){
    k=n-setpro;
    pb[0]=qb[0]=pb[1]=qb[1]=0.0;
    for (i=0;i<lena;i++){
     pb[i]=1.0/lena;
     for (j=0;j<lenb;j++){
        qb[j]=1.0/lenb;
        s=0.0;
        for (ii=2;ii<MAXE;ii++){
             s+= (PSICtemplate[0][i].frequency[ii] * PSICtarget[0][j].frequency[ii]);
        }
        score2[i][j] =  s ;
      }
    }
    lambda2[n]=lambda(pb,qb,score2,lena,lenb);
    printf("lambda2 %d = %f\n",n,lambda2[n]);
   }



   for (i=0;i<lena;i++){
   for (j=0;j<lenb;j++){
    matrix[j][i]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
      for (ii=0;ii<MAXE;ii++){
            matrix[j][i] +=  cte[n-1]*  (PSICtemplate[n][i].frequency[ii] * PSICtarget[n][j].frequency[ii]) * lambda2[n]/lambda1[n] ;
      }
     }
    }

    for (k=0;k<setmat;k++){
     for (ii=0;ii<MAXE;ii++){
            matrix[j][i] += ( 1.0 - ctp ) * weight[k] *   (PSICtemplate[0][i].frequency[ii] * PSICtarget[0][j].frequency[ii])* lambda2[n]/lambda1[n] ;
     }
    }
   }}

   for (i=0;i<lena;i++){
   for (j=0;j<lenb;j++){
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
      g[i]      +=  cte[n-1]*PSICtemplate[n][i].frequency[ii] * data[1][n-1][0][ii];
      ge[i]     +=  cte[n-1]*PSICtemplate[n][i].frequency[ii] * data[1][n-1][1][ii];
     }
     }
    }
    for (k=0;k<setmat;k++){
     for (ii=0;ii<MAXE;ii++){
      g[i]      += ( 1.0 - ctp ) * weight[k] * PSICtemplate[0][i].frequency[ii] * data[0][k][0][ii];
      ge[i]     += ( 1.0 - ctp ) * weight[k] * PSICtemplate[0][i].frequency[ii] * data[0][k][1][ii];
     }
    }
   }
   for (j=0;j<lenb;j++){
    p[j] =0.0;
    pe[j]=0.0;
    if ( checkpro==1 ){
     for (n=1;n<setpro;n++){
     for (jj=0;jj<MAXE;jj++){
      p[j]      +=  cte[n-1]*PSICtarget[n][j].frequency[jj] *data[1][n-1][0][jj];
      pe[j]     +=  cte[n-1]*PSICtarget[n][j].frequency[jj] *data[1][n-1][1][jj];
     }
     }
    }
    for (k=0;k<setmat;k++){
    for (jj=0;jj<MAXE;jj++){
      p[j]      += ( 1.0 - ctp ) * weight[k] * PSICtarget[0][j].frequency[jj] *data[0][k][0][jj];
      pe[j]     += ( 1.0 - ctp ) * weight[k] * PSICtarget[0][j].frequency[jj] *data[0][k][1][jj];
    }
    }
   }

   free_fvector(lambda1,0,setpro+setmat);
   free_fvector(lambda2,0,setpro+setmat);
   free_fvector(pb,0,MAXE+lena);
   free_fvector(qb,0,MAXE+lenb);
   free_Fmatrix(score1,0,MAXE,0,MAXE);
   free_Fmatrix(score2,0,lena,0,lenb);

   return  media;
 }
