#include  "inc.h"

int  EVDisland(m,a,b,g,ge,p,pe,lms,parameter,limit,tolerance)
float      **m;
sequence    a;
sequence    b;
float       *g,*ge,*p,*pe;
int         lms;
evd         *parameter;
float       limit;
float      *tolerance;
{
 float     **matrix;
 float       *gap,*gex,*pap,*pex;
 sequence    template,target;
 int         i,j,n;

 int         EVDlocal();
 float     **Fmatrix();
 float      *fvector();
 void        free_Fmatrix();
 void        free_vector();

  template.length=a.length+1;
  target.length  =b.length+1;
  strcpy(template.title,a.title);
  strcpy(target.title,b.title);


  matrix=Fmatrix(0,target.length,0,template.length);
  for (i=1;i<template.length;i++){
   for (j=1;j<target.length;j++){
      matrix[j][i]=m[j-1][i-1];
   }
  }
     gap     =fvector(0,template.length);
     gex     =fvector(0,template.length);
     pap     =fvector(0,target.length);
     pex     =fvector(0,target.length);


  for (i=1;i<template.length;i++){
      gap[i] = g[i];
      gex[i] = ge[i];
      matrix[0][i]=m[0][i-1];
      template.element[i]=a.element[i-1];
  }
  for (j=1;j<target.length;j++)  {
      pap[j] = p[j];
      pex[j] = pe[j];
      matrix[j][0]=m[j-1][0];
      target.element[j]=b.element[j-1];
  }

  matrix[0][0]=m[0][0];
  gap[0]=g[0];
  gex[0]=ge[0];
  pap[0]=p[0];
  pex[0]=pe[0];
  template.element[0]=0;
  target.element[0]  =0;

  (*parameter).lambda=1.0;
  (*parameter).K     =1.0;
  (*parameter).sigma =1.0; 
  (*parameter).Rc    =1;
  (*parameter).Sc    =1.0;
  (*parameter).min   =0.0;

  n=EVDlocal(matrix,template,target,gap,gex,pap,pex,lms,parameter,limit,tolerance);

  free_Fmatrix(matrix,0,target.length,0,template.length);
  free_fvector(gap,0,template.length);
  free_fvector(gex,0,template.length);
  free_fvector(pap,0,target.length);
  free_fvector(pex,0,target.length);

  return n;

}


int  EVDlocal(m,a,b,g,ge,p,pe,lms,parameter,limit,tolerance)
float      **m;
sequence    a;
sequence    b;
float       *g,*ge,*p,*pe;
int         lms;
evd        *parameter;
float       limit;
float      *tolerance;
{

/* work variables */

 int        i,j,ig,jg,up_e,left_e,l,ll,lk,found,locus,lsize;
 int        n,imin,locals;
 int        sumR;
 float      matchv,upv,leftv,data[6],score,matchs,ups,lefts;
 float     *xcore,cxcore,delta,max,min,minimum,cut_off;
 max3       maxima;
 evd        island[100];
 subalign   suba[DIM];
 

/* special variables */

 int          dimension,*dimensionj;   /* Number of possible alignments */
 tabulation **table;                  /* Table of edit distances */
 boundary    *border;                 /*Boundary limits of locals */

 
/* functions */

 float       *fvector();
 float        extension();
 float        ffmin();
 float      **Fmatrix();
 tabulation **Tmatrix();
 void         LocalTrace();
 max3         MaxData3();
 boundary    *bvector();

/* free functions */

 void         free_vector();
 void         free_bvector();
 void         free_Fmatrix();
 void         nrerror();

 
/*
*  Create Bottom-up table of edit and pointers
*/

 
 printf("Local Dynamic Algorithm Data to Calculate E-value parameters:\
 \n\t Sublocal Value:               %6.1e\
 \n\t Sublocal Size:                %7d\
 \n\t Scaling  Error:               %6.1e\
 \n\t Scaling  Factor:              %6.1e\
 \n\t Minimum  Score:               %6.1e\
 \n\t Level of Subalignment:        %6.1e\
 \n\t Limiting penalty extension:   %6.1e\
 \n\t Limiting gap extension:       %6.1e\n",limit,lms,tolerance[1],tolerance[0],tolerance[3],tolerance[4],tolerance[6],tolerance[2]);

 if (limit<=0.0){printf("ERROR: NULL Score Matrix\n"); dimension=1; return dimension;}

 table=Tmatrix(0,b.length,0,a.length);

 for(ig=0;ig<DIM;ig++){
   suba[ig].group=0;
   suba[ig].Imax=0;
   suba[ig].Jmax=0;
   for(jg=0;jg<MAXG;jg++){
    suba[ig].Ipair[jg]=0; 
    suba[ig].Jpair[jg]=0; 
   }
 }

  locals=0;                             /* Local alignments    */
  table[0][0].value=ffmin(m[0][0],limit);
  table[0][0].score=m[0][0];
  table[0][0].up=0;
  table[0][0].left=0;
  if (m[0][0]>0){table[0][0].match=1;
                 suba[locals].group=1;
                 suba[locals].Ipair[0]=0;
                 suba[locals].Jpair[0]=0;
                 locals++;}
  else          {table[0][0].match=0;}
  table[0][0].backu=0;
  table[0][0].backl=0;
  table[0][0].dimension=1;

  for (i=1;i<a.length; i++){ 
   table[0][i].value=ffmin(m[0][i],limit);
   table[0][i].score=m[0][i];
   table[0][i].up=0;
   table[0][i].left=0;
   if (m[0][i]>0){table[0][i].match=1;
                 suba[locals].group=1;
                 suba[locals].Ipair[0]=i;
                 suba[locals].Jpair[0]=0;
                 locals++;}
   else          {table[0][i].match=0;}
   table[0][i].backu=0;
   table[0][i].backl=0;
   table[0][i].dimension=1;
  }
  for (j=1;j<b.length; j++){
   table[j][0].value=ffmin(m[j][0],limit);
   table[j][0].score=m[j][0];
   table[j][0].up=0;
   table[j][0].left=0;
   if (m[j][0]>0){table[j][0].match=1;
                 suba[locals].group=1;
                 suba[locals].Ipair[0]=0;
                 suba[locals].Jpair[0]=j;
                 locals++;}
   else          {table[j][0].match=0;}
   table[j][0].backu=0;
   table[j][0].backl=0;
   table[j][0].dimension=1;
  }

 for (j=1;j<b.length; j++){
   for (i=1;i<a.length; i++){
         matchv  =table[j-1][i-1].value;
         matchs  =table[j-1][i-1].score;
         up_e    =table[j-1][i].backu;
         upv     =table[j-up_e-1][i].value;
         ups     =table[j-up_e-1][i].score;
         left_e  =table[j][i-1].backl;
         leftv   =table[j][i-left_e-1].value;
         lefts   =table[j][i-left_e-1].score;
         data[0] = ffmin((matchv + m[j][i]),limit);
         data[1] = ffmin((leftv  + m[j][i]),limit) - g[i] - extension(ge,i,left_e,tolerance);
         data[2] = ffmin((upv    + m[j][i]),limit) - p[j] - extension(pe,j,up_e,tolerance);
         data[3] = matchs + m[j][i];
         data[4] = lefts  + m[j][i] - g[i] - extension(ge,i,left_e,tolerance);
         data[5] = ups    + m[j][i] - p[j] - extension(pe,j,up_e,tolerance);
         maxima=MaxData3(data,limit,tolerance);
         table[j][i].value= maxima.value;
         table[j][i].score= maxima.score;
         table[j][i].match     = maxima.first;
         table[j][i].left      = maxima.second;
         table[j][i].up        = maxima.third;
         table[j][i].dimension =  table[j-1][i-1].dimension*maxima.first
                                + table[j][i-1].dimension*maxima.second
                                + table[j-1][i].dimension*maxima.third;
         if (table[j][i].dimension<=0){table[j][i].dimension=1;}
         if (  (maxima.first==1 && maxima.second==0 && maxima.third==0)
             ||(maxima.first==0 && maxima.second==0 && maxima.third==0)
            )
         {
            table[j][i].backu=0;
            table[j][i].backl=0;
         }else{
            table[j][i].backl = left_e + maxima.second;
            table[j][i].backu = up_e   + maxima.third;
         }
         if (
                (table[j][i].value > 0.0) &&
             (  (matchv > 0 && table[j][i].match==1)
              ||( leftv > 0 && table[j][i].left==1)
              ||(   upv > 0 && table[j][i].up==1) )
            ){
          if (locals>DIM){printf("too many local alignments (%d), exceed DIM (%d)\n",locals,DIM);locals=DIM-1;}
          if (locals>0){
            for (l=0;l<locals;l++){
             found=0;
             for (ll=0;ll < suba[l].group;ll++){
               if(
                      (    (i-1)==suba[l].Ipair[ll]
                       &&  (j-1)==suba[l].Jpair[ll]
                       &&  table[j][i].match==1
                      )
                   ||
                      (    (i-left_e-1)==suba[l].Ipair[ll]
                       &&  ( j )==suba[l].Jpair[ll]
                       &&  table[j][i].left==1
                      )
                   ||
                      (    ( i )==suba[l].Ipair[ll]
                       &&  (j-up_e-1)==suba[l].Jpair[ll]
                       &&  table[j][i].up==1
                      )
                 ){
                  if ( (i-suba[l].Ipair[0])>lms && (j-suba[l].Jpair[0])>lms){
                    score=  m[suba[l].Jpair[0]][suba[l].Ipair[0]]
                          + table[j][i].score
                          - table[suba[l].Jpair[0]][suba[l].Ipair[0]].score;
                    if (score>tolerance[3]){ found=1; locus=l; break; }
                  }else{found=1; locus=l; break; }
               }
             }
             if (found==1){break;}
            }
            if  (found==1){
             lsize=suba[locus].group;
             if (lsize < MAXG-1){
               suba[locus].group=lsize+1;
               suba[locus].Ipair[lsize]=i;
               suba[locus].Jpair[lsize]=j;
               if (suba[locus].Imax<i)  suba[locus].Imax=i;
               if (suba[locus].Jmax<j)  suba[locus].Jmax=j;
             }else{
               printf("Reached Maximum value of Long-Local Subalignment (MAXG)\n");
               if ( table[j][i].match==1 ) {
                  suba[locals].group=1;
                  suba[locals].Ipair[0]=i;
                  suba[locals].Jpair[0]=j;
                  suba[locals].Imax=i;
                  suba[locals].Jmax=j;
                  locals++;
               }
             }
            }else{
             suba[locals].group=1;
             suba[locals].Ipair[0]=i;
             suba[locals].Jpair[0]=j;
             suba[locals].Imax=i;
             suba[locals].Jmax=j;
             locals++;
            }
          }else{
            suba[locals].group=1;
            suba[locals].Ipair[0]=i;
            suba[locals].Jpair[0]=j;
            suba[locals].Imax=i;
            suba[locals].Jmax=j;
            locals++;
          }
         }
   }
 }

  if (locals>DIM){printf("too many local alignments (%d), exceed DIM (%d)\n",locals,DIM);locals=DIM-1;}



  sumR=0;
  xcore=fvector(0,DIM);
  border=bvector(0,1);
  max=tolerance[3];
  min=MAXL;
  for (l=0;l<locals;l++){
      lsize=suba[l].group-1;
      border[0].imax=suba[l].Imax;
      border[0].jmax=suba[l].Jmax;
      border[0].imin=suba[l].Ipair[0];
      border[0].jmin=suba[l].Jpair[0];
      if (border[0].jmin >=0 && border[0].imin >=0 && border[0].jmax<b.length && border[0].imax<a.length){
            score=   m[border[0].jmin][border[0].imin]
                   + table[border[0].jmax][border[0].imax].score
                   - table[border[0].jmin][border[0].imin].score;
            xcore[l]=score;
            if (score>max)max=score;
            if (score<min)min=score;
      }
      if  (score > tolerance[3] ) { sumR++; }
  }

  if (EVAL==1){
    minimum=MAXL;
    imin=0;
    delta=(max -min)/100;
    for (i=0;i<100;i++){
      island[i].lambda=1.0;
      island[i].K=1.0;
      island[i].sigma=1.0;
      cut_off =min + i*delta;
      cxcore=0.0;
      n=0;
      for (l=0;l<locals;l++) if ( xcore[l] >= cut_off) n++;
      for (l=0;l<locals;l++) if ( xcore[l] >= cut_off) cxcore += (xcore[l]-cut_off);
      if (n>0) cxcore= cxcore/n ;
      if (cxcore>0) island[i].lambda=log(1.0 + 1.0/cxcore);
      if (island[i].lambda<minimum) {minimum=island[i].lambda;imin=i;}
      if (n>0 && island[i].lambda>0 ) island[i].sigma=(exp(island[i].lambda) -1.0)/(island[i].lambda * sqrtf(exp(island[i].lambda)) * sqrtf((float)n) );
      island[i].K  = n * exp( cut_off * island[i].lambda) / (a.length*b.length) ;
      island[i].Rc = n;
      island[i].Sc = cxcore;
      island[i].min= min;
    }
    printf("Island method to calculate EVD parameters\n");
    printf("\t%s \t%s \t%s \t\t%s \t\t%s \t\t%s \n","Cut_off","Rc","Sc","Lambda","K","Sigma");
    for (i=0;i<100;i+=5){ 
      if (i>=imin && i<imin+5) printf("*");
      printf("\t%e\t%d\t%e\t%e\t%e\t%e\n",i*delta+min,island[i].Rc,island[i].Sc,island[i].lambda,island[i].K,island[i].sigma);
    }
    (*parameter)=island[imin];
  }
  (*parameter).min=min;

 free_bvector(border,0,1);
 free_fvector(xcore,0,DIM);
 free_Tmatrix(table,0,b.length,0,a.length);

 return  sumR;

}

