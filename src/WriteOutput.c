#include  "inc.h"

void WriteOutput(OUT,setpro,method,format,traduce,template,target,n,matrix,result)
FILE        *OUT;
int         setpro,method,format,n;
profile     *template,*target;
char        traduce[MAXE];
float       **matrix;
alignment   *result;
{

 int     fmt,i,j,kk,jj,ii,itmp0,itmp1,itgt0,itgt1,LenTmp,LenTgt;
 int     itmp0r,itmp1r,itgt0r,itgt1r,ir,nid,homo;
 int     **itmpir,**itgtir,count[MAXE],itmp,itgt,jtmp,jtgt,length;
 char    **seq_tmp,**seq_trg,*seq_cmp,**NameTmp,**NameTgt,*NameTmp0,*NameTgt0;
 char    **word_tmp,**word_trg,*word_cmp;
 float   mu,max,ident,homol;

 char    **cmatrix();
 char    *cvector();
 int     **imatrix();
 int     Consensus();
 void    free_cmatrix();
 void    free_cvector();
 void    free_imatrix();

  seq_tmp =cmatrix(0,template[0].size,0,MAXS);
  seq_trg =cmatrix(0,target[0].size,0,MAXS);
  seq_cmp =cvector(0,MAXS);
  word_tmp =cmatrix(0,template[0].size,0,MAXLINE);
  word_trg =cmatrix(0,target[0].size,0,MAXLINE);
  word_cmp =cvector(0,MAXLINE);

  itmpir =imatrix(0,template[0].size,0,MAXS);
  itgtir =imatrix(0,target[0].size,0,MAXS);

  NameTmp =cmatrix(0,template[0].size,0,MAXTITLE);
  NameTgt =cmatrix(0,target[0].size,0,MAXTITLE);
  NameTmp0=cvector(0,MAXTITLE);
  NameTgt0=cvector(0,MAXTITLE);


     memset(NameTmp0,'\0',MAXTITLE);
     memset(NameTgt0,'\0',MAXTITLE);
     for (jj=0;jj<template[0].size;jj++){ memset(NameTmp[jj],'\0',MAXTI);}
     for (jj=0;jj<target[0].size;jj++){ memset(NameTgt[jj],'\0',MAXTI);}

     LenTmp=strlen(template[0].main.title)<=MAXTITLE?strlen(template[0].main.title):MAXTITLE;
     LenTgt=strlen(target[0].main.title)<=MAXTITLE?strlen(target[0].main.title):MAXTITLE;
     strncpy(NameTmp0,template[0].main.title,LenTmp-1);
     strncpy(NameTgt0,target[0].main.title,LenTgt-1);

     if (format>9){fmt=format-10;}else{fmt=format;}
     mu =0.0;
     max=0.0;
     if (method==0){
      kk =0;
      for(jj=0;jj<target[0].main.length;jj++){for (ii=0;ii<template[0].main.length;ii++){
        max = max>=matrix[jj][ii]?max:matrix[jj][ii];
      }}
      for(jj=0;jj<target[0].main.length;jj++){for (ii=0;ii<template[0].main.length;ii++){
        if (matrix[jj][ii]>0.0 && matrix[jj][ii]<=max) {mu  += matrix[jj][ii]; kk++;}
      }}
      if (kk>0){mu=mu/kk;}
      mu=mu/3.0;
      max=max/5.0;
     }

     fprintf(OUT,"\n************************************\n");
     fprintf(OUT,"\nCOMPARE (Method %2d) ( %d x %d Sequences ) %s x %s  => %d alignments \n",
             method,template[0].size,target[0].size,NameTmp0,NameTgt0,n);

     for (i=0;i<n;i++)
     {

     for (kk=0;kk<setpro;kk++)
     {
      if (format>9) {fprintf(OUT,"\nPROPERTY %5d\n",kk);}
      else           {if (kk>0) {break;}                 }
 
      for (jj=0;jj<template[0].size;jj++){
       memset(NameTmp[jj],'\0',MAXTI);
       LenTmp=strlen(template[0].sequence[jj].title)<=MAXTI?strlen(template[0].sequence[jj].title):MAXTI;
       strncpy(NameTmp[jj],template[0].sequence[jj].title,LenTmp-1);
       strncat(NameTmp[jj],"               ",MAXTI-LenTmp);
       memset(seq_tmp[jj],'\0',MAXS);
      }

      for (jj=0;jj<target[kk].size;jj++){
       memset(NameTgt[jj],'\0',MAXTI);
       LenTgt=strlen(target[0].sequence[jj].title)<=MAXTI?strlen(target[0].sequence[jj].title):MAXTI;
       strncpy(NameTgt[jj],target[0].sequence[jj].title,LenTgt-1);
       strncat(NameTgt[jj],"               ",MAXTI-LenTgt);
       memset(seq_trg[jj],'\0',MAXS);
      }

      memset(seq_cmp,'\0',MAXS);
      nid=0;
      homo=0;
      for (j=0;j<result[i].length;j++)
      {

        itmp=result[i].template.element[j]-1;
        itgt=result[i].target.element[j]-1;

        itmp0 =result[i].Itemplate.position[j];
        itgt0 =result[i].Itarget.position[j];
        itmp0r=result[i].Itemplate.image[itmp0];
        itgt0r=result[i].Itarget.image[itgt0];
       

        jtmp=jtgt=0;

        for (jj=0;jj<MAXE;jj++){count[jj]=0;}
        for (jj=0;jj<template[kk].size;jj++){
            if (itmp>=0  ) {
             seq_tmp[jj][j]=traduce[template[kk].sequence[jj].element[itmp]];
             if( strncmp(&seq_tmp[jj][j],"-",1)){
              itmpir[jj][j]=template[kk].sequence[jj].position[itmp];
              count[template[kk].sequence[jj].element[itmp]]++;
             }else{
              if (j>0){itmpir[jj][j]=itmpir[jj][j-1];}else{itmpir[jj][j]=1;}
             }
            }else{
             seq_tmp[jj][j]='-';
             if (j>0){itmpir[jj][j]=itmpir[jj][j-1];}else{itmpir[jj][j]=1;}
            }
        }
        for (jj=0;jj<MAXE;jj++){
             count[jj]=100*count[jj]/template[kk].size;
             if (ConsensusID(count[jj],template[kk].size)==1){jtmp=jj;break;}
             //if (count[jj]>50){jtmp=jj;break;}
        }
            
        for (jj=0;jj<MAXE;jj++){count[jj]=0;}
        for (jj=0;jj<target[kk].size;jj++){
            if (itgt>=0){
             seq_trg[jj][j]=traduce[target[kk].sequence[jj].element[itgt]];
             if (strncmp(&seq_trg[jj][j],"-",1)){
              itgtir[jj][j]=target[kk].sequence[jj].position[itgt];
              count[target[kk].sequence[jj].element[itgt]]++;
             }else{
              if (j>0){itgtir[jj][j]=itgtir[jj][j-1];}else{itgtir[jj][j]=1;}
             }
            }else{
             seq_trg[jj][j]='-';
             if (j>0){itgtir[jj][j]=itgtir[jj][j-1];}else{itgtir[jj][j]=1;}
            }
        }
        for (jj=0;jj<MAXE;jj++){
             count[jj]=100*count[jj]/target[kk].size;
             if (ConsensusID(count[jj],target[kk].size)==1){jtgt=jj;break;}
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
      itmp0=result[i].Itemplate.position[0];
      itmp1=result[i].Itemplate.position[result[i].length-1];
      itgt0=result[i].Itarget.position[0];
      itgt1=result[i].Itarget.position[result[i].length-1];


      ir=0;
      while (itmp0 < 0 ) {ir++; itmp0=result[i].Itemplate.position[ir];}
      itmp0r=itmp0;
      ir=result[i].length-1;
      while (itmp1 < 0 ) {ir--; itmp1=result[i].Itemplate.position[ir];}
      itmp1r=itmp1;
      ir=0;
      while (itgt0 < 0 ) {ir++; itgt0=result[i].Itarget.position[ir];}
      itgt0r=itgt0;
      ir=result[i].length-1;
      while (itgt1 < 0 ) {ir--; itgt1=result[i].Itarget.position[ir];}
      itgt1r=itgt1;

      ident=100.0*nid/result[i].length;
      homol=100.0*(nid+homo)/result[i].length;

      fprintf(OUT,"\nMethod (%2d) Alignment (%d): %s_[%d-%d] versus %s_[%d-%d] (Score=%e  E-value=%e P-value=%e ID: %7.2f  Homology: %7.2f  )\n\n",
             method,i+1,NameTmp0,itmp0r,itmp1r,NameTgt0,itgt0r,itgt1r,
             result[i].score,result[i].evalue,result[i].pvalue,ident,homol);


      switch (fmt) {

        case 0:
          length=0;
          for (jj=0;jj<template[0].size;jj++){memset(word_tmp[jj],'\0',MAXLINE);}
          for (jj=0;jj<target[0].size;jj++){memset(word_trg[jj],'\0',MAXLINE);}
          memset(word_cmp,'\0',MAXLINE);
          for (j=0;j<result[i].length ;j++)
          {
            if (length==MAXLINE){
               fprintf(OUT,"\n");
               for (jj=0;jj<template[0].size;jj++){
                 fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTmp[jj],itmpir[jj][j-length],word_tmp[jj],itmpir[jj][j-1]);
               }
               fprintf(OUT,"                           %s         \n",word_cmp);
               for (jj=0;jj<target[0].size;jj++){
                 fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTgt[jj],itgtir[jj][j-length],word_trg[jj],itgtir[jj][j-1]);
               }
               length=0;
               for (jj=0;jj<template[0].size;jj++){memset(word_tmp[jj],'\0',MAXLINE);}
               for (jj=0;jj<target[0].size;jj++){memset(word_trg[jj],'\0',MAXLINE);}
               memset(word_cmp,'\0',MAXLINE);
            }      
            for (jj=0;jj<template[0].size;jj++){ strncat(word_tmp[jj],&seq_tmp[jj][j],1); }
            for (jj=0;jj<target[0].size;jj++){ strncat(word_trg[jj],&seq_trg[jj][j],1); }
            strncat(word_cmp,&seq_cmp[j],1);
            length++;
          } 
          fprintf(OUT,"\n");
          for (jj=0;jj<template[0].size;jj++){ fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTmp[jj],itmpir[jj][j-length],word_tmp[jj],itmpir[jj][j-1]);}
          fprintf(OUT,"                           %s         \n",word_cmp);
          for (jj=0;jj<target[0].size;jj++){ fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTgt[jj],itgtir[jj][j-length],word_trg[jj],itgtir[jj][j-1]);}
         break;

        case 1:

          for (jj=0;jj<template[0].size;jj++){
            length=0;
            fprintf(OUT,">P1;%s\n",NameTmp[jj]);
            fprintf(OUT,"%s\n",NameTmp[jj]);
            for (j=0;j<result[i].length;j++)
            {
              if (length==MAXLINE){fprintf(OUT,"\n");length=0;}
              fprintf(OUT,"%c",seq_tmp[jj][j]);
              length++;
            }
            fprintf(OUT,"\n");  
          }
          for (jj=0;jj<target[0].size;jj++){
            length=0;
            fprintf(OUT,">P1;%s\n",NameTgt[jj]);
            fprintf(OUT,"%s\n",NameTgt[jj]);
            for (j=0;j<result[i].length;j++)
            {
              if (length==MAXLINE){fprintf(OUT,"\n");length=0;}
              fprintf(OUT,"%c",seq_trg[jj][j]);
              length++;
            }  
            fprintf(OUT,"\n");
          }
         break;

        case 2:

          length=0;
          for (jj=0;jj<template[0].size;jj++){memset(word_tmp[jj],'\0',MAXLINE);}
          for (jj=0;jj<target[0].size;jj++){memset(word_trg[jj],'\0',MAXLINE);}
          for (j=0;j<result[i].length;j++)
          {
            if (length==MAXLINE){
               fprintf(OUT,"\n");
               for (jj=0;jj<template[0].size;jj++){
                 fprintf(OUT,"%s    %s  \n",NameTmp[jj],word_tmp[jj]);
               }
               for (jj=0;jj<target[0].size;jj++){
                 fprintf(OUT,"%s    %s  \n",NameTgt[jj],word_trg[jj]);
               }
               length=0;
               for (jj=0;jj<template[0].size;jj++){memset(word_tmp[jj],'\0',MAXLINE);}
               for (jj=0;jj<target[0].size;jj++){memset(word_trg[jj],'\0',MAXLINE);}
            }      
            for (jj=0;jj<template[0].size;jj++){ strncat(word_tmp[jj],&seq_tmp[jj][j],1); }
            for (jj=0;jj<target[0].size;jj++){ strncat(word_trg[jj],&seq_trg[jj][j],1); }
            length++;
          } 
          fprintf(OUT,"\n");
          for (jj=0;jj<template[0].size;jj++){ fprintf(OUT,"%s    %s  \n",NameTmp[jj],word_tmp[jj]);}
          for (jj=0;jj<target[0].size;jj++)  { fprintf(OUT,"%s    %s  \n",NameTgt[jj],word_trg[jj]);}
         break;

         default:

          length=0;
          for (jj=0;jj<template[0].size;jj++){memset(word_tmp[jj],'\0',MAXLINE);}
          for (jj=0;jj<target[0].size;jj++){memset(word_trg[jj],'\0',MAXLINE);}
          memset(word_cmp,'\0',MAXLINE);
          for (j=0;j<result[i].length ;j++)
          {
            if (length==MAXLINE){
               fprintf(OUT,"\n");
               for (jj=0;jj<template[0].size;jj++){
                 fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTmp[jj],itmpir[jj][j-length],word_tmp[jj],itmpir[jj][j-1]);
               }
               fprintf(OUT,"                           %s         \n",word_cmp);
               for (jj=0;jj<target[0].size;jj++){
                 fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTgt[jj],itgtir[jj][j-length],word_trg[jj],itgtir[jj][j-1]);
               }
               length=0;
               for (jj=0;jj<template[0].size;jj++){memset(word_tmp[jj],'\0',MAXLINE);}
               for (jj=0;jj<target[0].size;jj++){memset(word_trg[jj],'\0',MAXLINE);}
               memset(word_cmp,'\0',MAXLINE);
            }      
            for (jj=0;jj<template[0].size;jj++){ strncat(word_tmp[jj],&seq_tmp[jj][j],1); }
            for (jj=0;jj<target[0].size;jj++){ strncat(word_trg[jj],&seq_trg[jj][j],1); }
            strncat(word_cmp,&seq_cmp[j],1);
            length++;
          } 
          fprintf(OUT,"\n");
          for (jj=0;jj<template[0].size;jj++){ fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTmp[jj],itmpir[jj][j-length],word_tmp[jj],itmpir[jj][j-1]);}
          fprintf(OUT,"                           %s         \n",word_cmp);
          for (jj=0;jj<target[0].size;jj++){ fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTgt[jj],itgtir[jj][j-length],word_trg[jj],itgtir[jj][j-1]);}
         break;


      }
     }/* Change property */ 
     }


    free_cmatrix(NameTmp,0,template[0].size,0,MAXTI);
    free_cmatrix(NameTgt,0,target[0].size,0,MAXTI);
    free_cvector(NameTmp0,0,MAXTITLE);
    free_cvector(NameTgt0,0,MAXTITLE);

    free_cmatrix(seq_tmp,0,template[0].size,0,MAXS);
    free_cmatrix(seq_trg,0,target[0].size,0,MAXS);
    free_cmatrix(word_tmp,0,template[0].size,0,MAXLINE);
    free_cmatrix(word_trg,0,target[0].size,0,MAXLINE);
    free_imatrix(itmpir,0,template[0].size,0,MAXS);
    free_imatrix(itgtir,0,target[0].size,0,MAXS);
    free_cvector(seq_cmp,0,MAXS);
    free_cvector(word_cmp,0,MAXLINE);


}

