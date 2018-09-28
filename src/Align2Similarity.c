#include  "inc.h"

void Align2Similarity(a,template,target,matrix,traduce,simil)
alignment a;
profile   template;
profile   target;
float   **matrix;
char      traduce[MAXE];
float     simil[];
{

 int     i,j,kk,jj,ii,itmp0,itmp1,itgt0,itgt1;
 int     itmp0r,itmp1r,itgt0r,itgt1r,ir,nid,homo;
 int     **itmpir,**itgtir,count[MAXE],itmp,itgt,jtmp,jtgt;
 char    **seq_tmp,**seq_trg,*seq_cmp;
 float   mu,max;

 char    **cmatrix();
 char    *cvector();
 int     **imatrix();
 int     ConsensusID();

 void    free_cmatrix();
 void    free_cvector();
 void    free_imatrix();

  seq_tmp =cmatrix(0,template.size,0,MAXS);
  seq_trg =cmatrix(0,target.size,0,MAXS);
  seq_cmp =cvector(0,MAXS);

  itmpir =imatrix(0,template.size,0,MAXS);
  itgtir =imatrix(0,target.size,0,MAXS);

     mu =0.0;
     max=0.0;
     kk =0;
     for(jj=0;jj<target.main.length;jj++){for (ii=0;ii<template.main.length;ii++){
        max = max>=matrix[jj][ii]?max:matrix[jj][ii];
     }}
     for(jj=0;jj<target.main.length;jj++){for (ii=0;ii<template.main.length;ii++){
        if (matrix[jj][ii]>0.0 && matrix[jj][ii]<=max) {mu  += matrix[jj][ii]; kk++;}
     }}
     if (kk>0){mu=mu/kk;}
     mu=mu/3;
     max=max/5.0;
 

      nid=0;
      homo=0;
      for (j=0;j<a.length;j++)
      {

        itmp=a.template.element[j]-1;
        itgt=a.target.element[j]-1;

        itmp0 =a.Itemplate.position[j];
        itgt0 =a.Itarget.position[j];
        itmp0r=a.Itemplate.image[itmp0];
        itgt0r=a.Itarget.image[itgt0];
       

        jtmp=jtgt=0;

        for (jj=0;jj<MAXE;jj++){count[jj]=0;}
        for (jj=0;jj<template.size;jj++){
            if (itmp>=0  ) {
             seq_tmp[jj][j]=traduce[template.sequence[jj].element[itmp]];
             if( strncmp(&seq_tmp[jj][j],"-",1)){
              itmpir[jj][j]=template.sequence[jj].position[itmp];
              count[template.sequence[jj].element[itmp]]++;
             }else{
              if (j>0){itmpir[jj][j]=itmpir[jj][j-1];}else{itmpir[jj][j]=1;}
             }
            }else{
             seq_tmp[jj][j]='-';
             if (j>0){itmpir[jj][j]=itmpir[jj][j-1];}else{itmpir[jj][j]=1;}
            }
        }
        for (jj=0;jj<MAXE;jj++){
             count[jj]=100*count[jj]/template.size;
             if (ConsensusID(count[jj],template.size)==1){jtmp=jj;break;}
             //if (count[jj]>50){jtmp=jj;break;}
        }
            
        for (jj=0;jj<MAXE;jj++){count[jj]=0;}
        for (jj=0;jj<target.size;jj++){
            if (itgt>=0){
             seq_trg[jj][j]=traduce[target.sequence[jj].element[itgt]];
             if (strncmp(&seq_trg[jj][j],"-",1)){
              itgtir[jj][j]=target.sequence[jj].position[itgt];
              count[target.sequence[jj].element[itgt]]++;
             }else{
              if (j>0){itgtir[jj][j]=itgtir[jj][j-1];}else{itgtir[jj][j]=1;}
             }
            }else{
             seq_trg[jj][j]='-';
             if (j>0){itgtir[jj][j]=itgtir[jj][j-1];}else{itgtir[jj][j]=1;}
            }
        }
        for (jj=0;jj<MAXE;jj++){
             count[jj]=100*count[jj]/target.size;
             if (ConsensusID(count[jj],target.size)==1){jtgt=jj;break;}
             //if (count[jj]>50){jtgt=jj;break;}
        }
            
        if (itgt0r >= 0 && itmp0r >= 0 
            && strncmp(&traduce[jtgt],"-",1)
            && strncmp(&traduce[jtmp],"-",1)){ 
          if (jtmp==jtgt){seq_cmp[j]=traduce[jtmp];nid++;}
          else{ if (   matrix[itmp0r][itgt0r] > mu 
                    && matrix[itmp0r][itgt0r] > max){seq_cmp[j]='+';homo++;}
                else                                {seq_cmp[j]=' ';}}}
        else                                        {seq_cmp[j]=' ';}

      }
/*
      itmp0=a.Itemplate.position[0];
      itmp1=a.Itemplate.position[a.length-1];
      itgt0=a.Itarget.position[0];
      itgt1=a.Itarget.position[a.length-1];


      ir=0;
      while (itmp0 < 0 ) {ir++; itmp0=a.Itemplate.position[ir];}
      itmp0r=itmp0;
      ir=a.length-1;
      while (itmp1 < 0 ) {ir--; itmp1=a.Itemplate.position[ir];}
      itmp1r=itmp1;
      ir=0;
      while (itgt0 < 0 ) {ir++; itgt0=a.Itarget.position[ir];}
      itgt0r=itgt0;
      ir=a.length-1;
      while (itgt1 < 0 ) {ir--; itgt1=a.Itarget.position[ir];}
      itgt1r=itgt1;
*/


      if (a.length>0)simil[0]=100.0*nid/a.length;
      if (a.length>0)simil[1]=100.0*(nid+homo)/a.length;
 
    free_cmatrix(seq_tmp,0,template.size,0,MAXS);
    free_cmatrix(seq_trg,0,target.size,0,MAXS);
    free_cvector(seq_cmp,0,MAXS);
    free_imatrix(itmpir,0,template.size,0,MAXS);
    free_imatrix(itgtir,0,target.size,0,MAXS);


}

