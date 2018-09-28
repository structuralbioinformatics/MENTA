#include  "inc.h"
/* Subroutine developed by B. Oliva  for MEnTA. */

int  GlobalDyna(m,a,b,g,ge,p,pe,evd_param,limit,tolerance,pssm,r)
float      **m;
sequence    a;
sequence    b;
float       *g,*ge,*p,*pe;
evd         evd_param;
float       limit;
float       tolerance[];
float     **pssm;
alignment  *r;
{
 float     **matrix,maxscc,minsc;
 float       *gap,*gex,*pap,*pex;
 sequence    template,target;
 int         i,j,n;
 int         GlobalDynamo();
 float     **Fmatrix();
 float      *fvector();
 void        free_Fmatrix();
 void        free_fvector();

  maxscc=0.0;
  minsc=MAXSC;

  template.length=a.length+1;
  target.length  =b.length+1;
  strcpy(template.title,a.title);
  strcpy(target.title,b.title);

  matrix  =Fmatrix(0,target.length,0,template.length);
  gap     =fvector(0,template.length);
  gex     =fvector(0,template.length);
  pap     =fvector(0,target.length);
  pex     =fvector(0,target.length);


  for (i=1;i<template.length;i++){
   for (j=1;j<target.length;j++){
      matrix[j][i]=m[j-1][i-1];
      if (matrix[j][i]>=maxscc){maxscc=matrix[j][i];}
      if (maxscc>MAXSC){printf("\n%CAUTION MAXSC < %f \n",maxscc);}
      if (matrix[j][i]<=minsc){minsc=matrix[j][i];}
   }
  }

 
  for (i=1;i<template.length;i++){
      matrix[0][i]=0.0  ;
      gap[i] = g[i-1];
      gex[i] = ge[i-1];
      template.element[i]=a.element[i-1];
  }
  for (j=1;j<target.length;j++)  {
      matrix[j][0]=0.0  ;
      pap[j] = p[j-1];
      pex[j] = pe[j-1];
      target.element[j]=b.element[j-1];
  }

  matrix[0][0]=0.0;
  matrix[0][1]=0.0;
  matrix[1][0]=0.0;

  gap[0]=g[0];
  gex[0]=ge[0];
  pap[0]=p[0];
  pex[0]=pe[0];
  template.element[0]=0;
  target.element[0]  =0;


  n=GlobalDynamo(matrix,template,target,gap,gex,pap,pex,evd_param,limit,tolerance,pssm,r);
  free_Fmatrix(matrix,0,target.length,0,template.length);
  free_fvector(gap,0,template.length);
  free_fvector(gex,0,template.length);
  free_fvector(pap,0,target.length);
  free_fvector(pex,0,target.length);
  return n;

}

int  GlobalDynamo(m,a,b,g,ge,p,pe,evd_param,limit,tolerance,pssm,r)
float      **m;
sequence    a;
sequence    b;
float       *g,*ge,*p,*pe;
evd         evd_param;
float       limit;
float       tolerance[];
float     **pssm;
alignment  *r;
{

/* work variables */

 int      i,j;
 int      up_e,left_e,u3d,l3d;
 float    size,matchv,upv,leftv,data[6],score,matchs,ups,lefts;
 max3     maximaUL;

/* special variables */

 int          dimension;   /* Number of possible alignments */
 tabulation **table;       /* Table of edit distances */

 int   dim_data,idata_matrix[5][MAXS][MAXS],data_seq[2][MAXS];
 float data_matrix[5][MAXS][MAXS];



 
/* functions */

 tabulation **Tmatrix();  
 alignment   *avector(); 
 int         *ivector();
 max3         MaxDataG3();
 void         free_avector();
 void         free_ivector();
 int          GlobalTrace();
 float        extension(),ffmin();
 double       evalue(),pvalue();
 void         ashell();

 void         PrintIMatrix();
 void         PrintFMatrix();



/*
*  Create Bottom-up table of edit and pointers
*/

 printf("Global Dynamic Algorithm Data:\
 \n\t Maximum Local Value:          %6.1e\
 \n\t Scaling  Error:               %6.1e\
 \n\t Scaling  Factor:              %6.1e\
 \n\t Minimum  Score:               %6.1e\
 \n\t Level of Subalignment:        %6.1e\
 \n\t Limiting penalty extension:   %6.1e\
 \n\t Limiting gap extension:       %6.1e\n",limit,tolerance[1],tolerance[0],tolerance[3],tolerance[4],tolerance[6],tolerance[2]);

 if (limit<=0.0){printf("ERROR: NULL Score Matrix\n"); dimension=0; return dimension;}

 table=Tmatrix(0,b.length,0,a.length);

  table[0][0].value=m[0][0];
  table[0][0].score=m[0][0];
  table[0][0].up=0;
  table[0][0].left=0;
  table[0][0].match=1;
  table[0][0].backu=0;
  table[0][0].backl=0;
  table[0][0].dimension=1;

  for (i=1;i<a.length; i++){ 
   table[0][i].value=m[0][i];
   table[0][i].score=m[0][i];
   table[0][i].up=0;
   table[0][i].left=1;
   table[0][i].match=0;
   table[0][i].backu=0;
   table[0][i].backl=i;
   table[0][i].dimension=1;
  }
  for (j=1;j<b.length; j++){
   table[j][0].value=m[j][0];
   table[j][0].score=m[j][0];
   table[j][0].up=1;
   table[j][0].left=0;
   table[j][0].match=0;
   table[j][0].backu=j;
   table[j][0].backl=0;
   table[j][0].dimension=1;
  }

 for (j=1;j<b.length; j++){
   for (i=1;i<a.length; i++){
         matchv  =table[j-1][i-1].value;
         matchs  =table[j-1][i-1].score;
     u3d=up_e    =table[j-1][i].backu;
         upv     =table[j-up_e-1][i].value;
         ups     =table[j-up_e-1][i].score;
     l3d=left_e  =table[j][i-1].backl;
         leftv   =table[j][i-left_e-1].value;
         lefts   =table[j][i-left_e-1].score;
         data[3] = matchs + m[j][i];
         data[4] = lefts  + m[j][i] - g[i] - extension(ge,i,left_e,tolerance);
         data[5] = ups    + m[j][i] - p[j] - extension(pe,j,up_e,tolerance);

 	 if (i==1||j==1 ){
         data[0] =  matchv + m[j][i];
         data[1] =  leftv  + m[j][i];
         data[2] =  upv    + m[j][i];
         data[3] =  matchs + m[j][i];
         data[4] =  lefts  + m[j][i];
         data[5] =  ups    + m[j][i];
	 }


         data[0] =ffmin((matchv + m[j][i]),limit);
         if (l3d>0){
             if (i-l3d>1){data[1] =ffmin((table[j][i-l3d-1].value + m[j][i]),limit) -g[i] -extension(ge,i,l3d,tolerance);}
             else        {data[1] =ffmin((table[j][i-l3d-1].value + m[j][i]),limit);}}
         else      {
             if (i>1)    {data[1] =ffmin((table[j][i-1].value + m[j][i]),limit) -g[i] ;}
             else        {data[1] =ffmin((table[j][i-1].value + m[j][i]),limit);}}
         if (u3d>0){
             if (j-u3d>1){data[2] =ffmin((table[j-u3d-1][i].value + m[j][i]),limit) -p[j] -extension(pe,j,u3d,tolerance);}
             else        {data[2] =ffmin((table[j-u3d-1][i].value + m[j][i]),limit);}}
         else      {
             if (j>1)    {data[2] =ffmin((table[j-1][i].value + m[j][i]),limit) -p[j];}
             else        {data[2] =ffmin((table[j-1][i].value + m[j][i]),limit);}}


         if (j==b.length-1){data[1] =ffmin(table[j][i-l3d-1].value,limit);}
         if (i==a.length-1){data[2] =ffmin(table[j-u3d-1][i].value,limit);}

         if (i==1||j==1 ){
          data[1] = ffmin(table[j][i-1].value,limit);
          data[2] = ffmin(table[j-1][i].value,limit);
	 }

         maximaUL=MaxDataG3(data,limit,tolerance);
         table[j][i].value = maximaUL.value;
         table[j][i].score = maximaUL.score;

         table[j][i].match     = maximaUL.first;
         table[j][i].left      = maximaUL.second;
         table[j][i].up        = maximaUL.third;
         table[j][i].dimension =  table[j-1][i-1].dimension*maximaUL.first 
                                + table[j][i-1].dimension*maximaUL.second
                                + table[j-1][i].dimension*maximaUL.third;
         if (table[j][i].dimension > DIM ) {
            printf("allocation problem when building up Table in GLOBALDYNA,DIMENSION[%d,%d]= %d \n",j,i,table[j][i].dimension);
         }
         if (maximaUL.first==1 && maximaUL.second==0 && maximaUL.third==0 ) {
            table[j][i].backu=0;
            table[j][i].backl=0;
         }else{
            table[j][i].backl = l3d + maximaUL.second;
            table[j][i].backu = u3d + maximaUL.third;
         }
   }
 }

 
 dimension=table[b.length-1][a.length-1].dimension;
 score    =table[b.length-1][a.length-1].score;


/*****************

      for (j=0;j<b.length; j++){
        data_seq[0][j]=b.element[j];
      for (i=0;i<a.length; i++){
        data_seq[1][i]=a.element[i];
        dim_data=3;
        data_matrix[0][j][i]=m[j][i];
        data_matrix[1][j][i]=table[j][i].value;
        data_matrix[2][j][i]=table[j][i].score;
        idata_matrix[0][j][i]=table[j][i].match;
        idata_matrix[1][j][i]=table[j][i].up;
        idata_matrix[2][j][i]=table[j][i].left;
      }}
      PrintFMatrix(a.length,b.length,dim_data,data_matrix,data_seq);
      PrintIMatrix(a.length,b.length,dim_data,idata_matrix,data_seq);

*****************/
 size= a.length * b.length;
 if (dimension>DIM){ printf("allocation problem in GLOBALDYNA, DIMENSION= %d \n",dimension);dimension=DIM-1; }
 if ( dimension < DIM){
/*Traceback*/  
   dimension=GlobalTrace(m,a,b,g,ge,p,pe,tolerance,table,dimension,pssm,r);
   for (i=0;i<dimension;i++){
            r[i].evalue=evalue(evd_param,r[i],size);
            r[i].pvalue=pvalue(evd_param,r[i],size);
   }
   ashell(dimension,r);
 }else{dimension=0;}

 free_Tmatrix(table,0,b.length,0,a.length);

 return dimension;

}

