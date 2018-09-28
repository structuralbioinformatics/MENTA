#include  "inc.h"

void WriteProfile(OUT,traduce,template)
FILE        *OUT;
profile     *template;
char        traduce[MAXE];
{

 int     i,j,jj,itmp0,itmp1,LenTmp;
 int     itmp0r,itmp1r,ir,nid,homo;
 int     **itmpir,count[MAXE],itmp,jtmp,length;
 char    **seq_tmp,*seq_cmp,**NameTmp,*NameTmp0;
 char    **word_tmp,*word_cmp;

 char    **cmatrix();
 char    *cvector();
 int     **imatrix();
 void    free_cmatrix();
 void    free_cvector();
 void    free_imatrix();

  seq_tmp =cmatrix(0,template[0].size,0,MAXS);
  seq_cmp =cvector(0,MAXS);
  word_tmp =cmatrix(0,template[0].size,0,MAXLINE);
  word_cmp =cvector(0,MAXLINE);

  itmpir =imatrix(0,template[0].size,0,MAXS);

  NameTmp =cmatrix(0,template[0].size,0,MAXTI);
  NameTmp0=cvector(0,MAXTI);


     memset(NameTmp0,'\0',MAXTI);
     for (jj=0;jj<template[0].size;jj++){ memset(NameTmp[jj],'\0',MAXTI);}

     LenTmp=strlen(template[0].main.title)<=MAXTI?strlen(template[0].main.title):MAXTI;
     strncpy(NameTmp0,template[0].main.title,LenTmp-1);
     fprintf(OUT,"\n******************************\n");
     fprintf(OUT,"\n# PROFILE %s  \n", template[0].main.title);


      for (jj=0;jj<template[0].size;jj++){
       memset(NameTmp[jj],'\0',MAXTI);
       LenTmp=strlen(template[0].sequence[jj].title)<=MAXTI?strlen(template[0].sequence[jj].title):MAXTI;
       strncpy(NameTmp[jj],template[0].sequence[jj].title,LenTmp-1);
       strncat(NameTmp[jj],"               ",MAXTI-LenTmp);
       memset(seq_tmp[jj],'\0',MAXS);
      }
      memset(seq_cmp,'\0',MAXS);
      nid=0;
      for (j=0;j<template[0].main.length;j++)
      {
        itmp=template[0].main.element[j]-1;
        for (jj=0;jj<MAXE;jj++){count[jj]=0;}
        for (jj=0;jj<template[0].size;jj++){
            seq_tmp[jj][j]=traduce[template[0].sequence[jj].element[itmp]];
            itmpir[jj][j]=template[0].sequence[jj].position[itmp];
            count[template[0].sequence[jj].element[itmp]]++;
        }
        jtmp=1000;
        for (jj=0;jj<MAXE;jj++){
             count[jj]=100*count[jj]/template[0].size;
             if (count[jj]>50){jtmp=jj;break;}
        }
        if (jtmp < 1000){
         if (count[jtmp]<=50){strncpy(&seq_cmp[j],".",1);}else{seq_cmp[j]=traduce[jtmp];nid++;}
        }else{
         strncpy(&seq_cmp[j],".",1);
        }
      }
      length=0;
      for (jj=0;jj<template[0].size;jj++){memset(word_tmp[jj],'\0',MAXLINE);}
      memset(word_cmp,'\0',MAXLINE);
      for (j=0;j<template[0].main.length;j++)
      {
            if (length==MAXLINE){
               fprintf(OUT,"\n");
               for (jj=0;jj<template[0].size;jj++){
                 fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTmp[jj],itmpir[jj][j-length],word_tmp[jj],itmpir[jj][j-1]);
               }
               fprintf(OUT,"CONSENSUS_____             %s         \n",word_cmp);
               length=0;
               for (jj=0;jj<template[0].size;jj++){memset(word_tmp[jj],'\0',MAXLINE);}
               memset(word_cmp,'\0',MAXLINE);
            }      
            for (jj=0;jj<template[0].size;jj++){ strncat(word_tmp[jj],&seq_tmp[jj][j],1); }
            strncat(word_cmp,&seq_cmp[j],1);
            length++;
      } 
      fprintf(OUT,"\n");
      for (jj=0;jj<template[0].size;jj++){ fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTmp[jj],itmpir[jj][j-length],word_tmp[jj],itmpir[jj][j-1]);}
      fprintf(OUT,"CONSENSUS_____             %s         \n",word_cmp);
      
 


    free_cmatrix(NameTmp,0,template[0].size,0,MAXTI);
    free_cvector(NameTmp0,0,MAXTI);

    free_cmatrix(seq_tmp,0,template[0].size,0,MAXS);
    free_cmatrix(word_tmp,0,template[0].size,0,MAXLINE);
    free_cvector(word_cmp,0,template[0].size,0,MAXLINE);
    free_imatrix(itmpir,0,template[0].size,0,MAXS);
    free_cvector(seq_cmp,0,MAXS);


}

