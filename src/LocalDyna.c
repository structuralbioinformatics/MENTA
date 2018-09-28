#include  "inc.h"

int  LocalDyna(m,a,b,g,ge,p,pe,lms,evd_param,limit,tolerance,pssm,r)
float      **m;
sequence    a;
sequence    b;
float       *g,*ge,*p,*pe;
int         lms;
evd         evd_param;
float       limit;
float      *tolerance;
float     **pssm;
alignment  *r;
{
 float     **matrix;
 float       *gap,*gex,*pap,*pex;
 sequence    template,target;
 int         i,j,n;
 int         LocalDynamo();
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
  n=LocalDynamo(matrix,template,target,gap,gex,pap,pex,lms,evd_param,limit,tolerance,pssm,r);
  free_Fmatrix(matrix,0,target.length,0,template.length);
  free_fvector(gap,0,template.length);
  free_fvector(gex,0,template.length);
  free_fvector(pap,0,target.length);
  free_fvector(pex,0,target.length);

  return n;

}
/*  PROGRAM  MENTA:  Multiple Entities Alignment */
/*  Baldomero Oliva. Computational Structural Biology Laboratory */
/*  Universitat Pompeu Fabra */
/*  Barcelona. Catalonia, Spain (EU */


int  LocalDynamo(m,a,b,g,ge,p,pe,lms,evd_param,limit,tolerance,pssm,r)
float      **m;
sequence    a;
sequence    b;
float       *g,*ge,*p,*pe;
int         lms;
evd         evd_param;
float       limit;
float      *tolerance;
float     **pssm;
alignment  *r;
{

/* work variables */

 int        i,j,ii,jj,up_e,left_e,l,ll,lk,found,locus,lsize;
 int        locals;
 int        sumD;
 float      size,matchv,upv,leftv,data[6],score,dscore,dscore2,matchs,ups,lefts;
 float    **pssm_tmp;
 max3       maxima;
 alignment *rl,*avector();
 subalign   suba[DIM];


/* special variables */

 int          dimension,*dimensionj;   /* Number of possible alignments */
 tabulation **table;                  /* Table of edit distances */
 boundary    *border;                 /*Boundary limits of locals */

 
/* functions */

 double       evalue();
 double       pvalue();
 float       *fvector();
 float        extension();
 float        ffmin();
 float      **Fmatrix();
 tabulation **Tmatrix();
 void         LocalTrace();
 max3         MaxData3();
 boundary    *bvector();
 void        ashell();

/* free functions */

 void         free_vector();
 void         free_avector();
 void         free_bvector();
 void         free_Fmatrix();
 void         nrerror();

 int     itmp0r,itmp1r,itgt0r,itgt1r;
 int     itmp0,itmp1,itgt0,itgt1,ir;
 int     ig,jg;

/*
*  Create Bottom-up table of edit and pointers
*/


 printf("Local Dynamic Algorithm Data:\
 \n\t Sublocal Value:              %6.1e\
 \n\t Sublocal Size:               %7d\
 \n\t Scaling  Error:              %6.1e\
 \n\t Scaling  Factor:             %6.1e\
 \n\t Minimum  Score:              %6.1e\
 \n\t Level of Subalignment:       %6.1e\
 \n\t Limiting penalty extension:  %6.1e\
 \n\t Limiting gap extension:      %6.1e\n",limit,lms,tolerance[1],tolerance[0],tolerance[3],tolerance[4],tolerance[6],tolerance[2]);

 if (limit<=0.0){printf("ERROR: NULL Score Matrix\n"); dimension=0; return dimension;}

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
   if (locals>DIM){printf("STOP local search, too many local alignments (%d), exceed DIM (%d)\n",locals,DIM);break;}
   for (i=1;i<a.length; i++){
         if (locals>DIM){printf("STOP local search, too many local alignments (%d), exceed DIM (%d)\n",locals,DIM);break;}
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
                  if ( (i-suba[l].Ipair[0]) >lms && (j-suba[l].Jpair[0])>lms){
                    score=  m[suba[l].Jpair[0]][suba[l].Ipair[0]]
                          + table[j][i].score
                          - table[suba[l].Jpair[0]][suba[l].Ipair[0]].score;
                    if (score>tolerance[3]){ found=1; locus=l; break; }
                    if (tolerance[3]<-1.0e+50){ found=1; locus=l; break; }
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


  pssm_tmp=Fmatrix(0,b.length+1,0,a.length+1);

  size=a.length * b.length;
  sumD=0;
  dscore=dscore2=0.0;
  border=bvector(0,1);
  rl=avector(0,DIM);
  for (l=0;l<locals;l++){
      lsize=suba[l].group-1;
      border[0].imax=suba[l].Imax;
      border[0].jmax=suba[l].Jmax;
      border[0].imin=suba[l].Ipair[0];
      border[0].jmin=suba[l].Jpair[0];
      if (sumD>=DIM){ printf("allocation failure in DYNA: INCREASE DIM \n"); sumD=DIM ; break;  }
      if  (  ( border[0].imax - border[0].imin ) > lms
           &&( border[0].jmax - border[0].jmin ) > lms )
      {
       dimension=1+(table[border[0].jmax][border[0].imax].dimension-table[border[0].jmin][border[0].imin].dimension);
       if (dimension>DIM){printf("Re-dimension Local Subalignment[%d] from %d to DIM \n",l,dimension);dimension=DIM-1;} 
       if (dimension > 0 &&  dimension<DIM) {
         if (sumD>=DIM){ printf("allocation failure in DYNA: INCREASE DIM \n"); break;  }
         score=   m[border[0].jmin][border[0].imin]
                + table[border[0].jmax][border[0].imax].score
                - table[border[0].jmin][border[0].imin].score;
         if ( (dimension<DIM && score > tolerance[3])  || (dimension<DIM && tolerance[3]<-1.0e+50) ){
          dimensionj=&dimension;
          for  (jj=0;jj<=b.length; jj++){ for (ii=0;ii<=a.length; ii++){ pssm_tmp[jj][ii]=0.0;}}
          LocalTrace(m,a,b,g,ge,p,pe,table,border,dimension,dimensionj,lms,tolerance,pssm_tmp,rl);
          if (*dimensionj > 0 ){
           for  (jj=1;jj<=b.length; jj++){ for (ii=1;ii<=a.length; ii++){ pssm[jj][ii]=pssm[jj][ii]+pssm_tmp[jj][ii]/(*dimensionj);}}
          }
          for (lk=0;lk<(*dimensionj);lk++){
           rl[lk].evalue= evalue(evd_param,rl[lk],size);
           rl[lk].pvalue= pvalue(evd_param,rl[lk],size);
           itmp0=rl[lk].Itemplate.position[0];
           itmp1=rl[lk].Itemplate.position[rl[lk].length-1];
           itgt0=rl[lk].Itarget.position[0];
           itgt1=rl[lk].Itarget.position[rl[lk].length-1];
           itmp0r=rl[lk].Itemplate.image[itmp0];
           itgt0r=rl[lk].Itarget.image[itgt0];
           itmp1r=rl[lk].Itemplate.image[itmp1];
           itgt1r=rl[lk].Itarget.image[itgt1];
           ir=0;
           while (itmp0 < 0 ) {ir++; itmp0=rl[lk].Itemplate.position[ir];}
           itmp0r=itmp0;
           ir=rl[lk].length-1;
           while (itmp1 < 0 ) {ir--; itmp1=rl[lk].Itemplate.position[ir];}
           itmp1r=itmp1;
           ir=0;
           while (itgt0 < 0 ) {ir++; itgt0=rl[lk].Itarget.position[ir];}
           itgt0r=itgt0;
           ir=rl[lk].length-1;
           while (itgt1 < 0 ) {ir--; itgt1=rl[lk].Itarget.position[ir];}
           itgt1r=itgt1;
           if ( (ir>0) && ((itmp1r-itmp0r + 1 )>lms) && ((itgt1r-itgt0r + 1 )>lms) && sumD < DIM ){
              r[sumD]=rl[lk];
              dscore  += score-rl[lk].score;
              dscore2 += dscore*dscore;
              sumD++;
           }
          }
         }
       }
      }
  }
  if (sumD>0) dscore =dscore/sumD;
  if (sumD>0) dscore2=dscore2/sumD;
  printf("Media-error on Score approach = %e | rmsd = %e\n",dscore,sqrtf(fabs(dscore2 - dscore*dscore)));
  ashell(sumD,r);

 free_avector(rl,0,DIM);
 free_bvector(border,0,1);
 free_Tmatrix(table,0,b.length,0,a.length);
 free_Fmatrix(pssm_tmp,0,b.length+1,0,a.length+1);

 return  sumD;

}

