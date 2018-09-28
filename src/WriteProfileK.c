#include  "inc.h"

void WriteProfileK(OUT,traduce,template,k,kk)
FILE        *OUT;
profile     **template;
char        traduce[MAXE];
int         k,kk;
{

 int     i,j,jj,itmp0,itmp1,LenTmp;
 int     itmp0r,itmp1r,ir,nid,homo;
 int     **itmpir,count[MAXE],itmp,jtmp,length;
 char    **seq_tmp,*seq_cmp,**NameTmp,*NameTmp0;
 char    **word_tmp,*word_cmp;

 char    **cmatrix();
 char    *cvector();
 int     **imatrix();
 int     ConsensusID();
 void    free_cmatrix();
 void    free_cvector();
 void    free_imatrix();


  seq_tmp =cmatrix(0,template[kk][k].size,0,MAXS);
  seq_cmp =cvector(0,MAXS);
  word_tmp =cmatrix(0,template[kk][k].size,0,MAXLINE);
  word_cmp =cvector(0,MAXLINE);

  itmpir =imatrix(0,template[kk][k].size,0,MAXS);

  NameTmp =cmatrix(0,template[kk][k].size,0,MAXTI);
  NameTmp0=cvector(0,MAXTI);


     memset(NameTmp0,'\0',MAXTI);
     for (jj=0;jj<template[kk][k].size;jj++){ memset(NameTmp[jj],'\0',MAXTI);}

     LenTmp=strlen(template[kk][k].main.title)<=MAXTI?strlen(template[kk][k].main.title):MAXTI;
     strncpy(NameTmp0,template[kk][k].main.title,LenTmp-1);
     fprintf(OUT,"\n******************************\n");
     fprintf(OUT,"\n# PROFILE [%d] %s  \n# PROFILE [%d] Include Cluster of Profiles: ", kk,template[kk][k].main.title,kk);
     for (jj=0;jj<template[kk][k].cluster_dimension;jj++){fprintf(OUT," %d",template[kk][k].cluster[jj]);}
     fprintf(OUT,"\n# PROFILE [%d] Score=%e  E-value=%e P-value=%e ID: %7.2f  Homology: %7.2f \n",kk,template[kk][k].score,template[kk][k].evalue,template[kk][k].pvalue,template[kk][k].ident,template[kk][k].homol);
     //fprintf(OUT,"# PROFILE [%d] Action %d\n",kk,template[kk][k].action);
     fprintf(OUT,"# PROFILE [%d] Number of Sequences %d\n",kk,template[kk][k].size);
     printf("\n******************************\n");
     printf("\n# PROFILE [%d] %s  \n# PROFILE [%d] Include Cluster of Profiles: ", kk,template[kk][k].main.title,kk);
     for (jj=0;jj<template[kk][k].cluster_dimension;jj++){printf(" %d",template[kk][k].cluster[jj]);}
     printf("\n# PROFILE [%d] Score=%e  E-value=%e P-value=%e ID: %7.2f  Homology: %7.2f \n",kk,template[kk][k].score,template[kk][k].evalue,template[kk][k].pvalue,template[kk][k].ident,template[kk][k].homol);
     //printf("# PROFILE [%d] Action %d\n",kk,template[kk][k].action);
     printf("# PROFILE [%d] Number of Sequences %d\n",kk,template[kk][k].size);


      for (jj=0;jj<template[kk][k].size;jj++){
       memset(NameTmp[jj],'\0',MAXTI);
       LenTmp=strlen(template[kk][0].sequence[jj].title)<=MAXTI?strlen(template[kk][0].sequence[jj].title):MAXTI;
       strncpy(NameTmp[jj],template[kk][0].sequence[jj].title,LenTmp-1);
       strncat(NameTmp[jj],"               ",MAXTI-LenTmp);
       memset(seq_tmp[jj],'\0',MAXS);
      }
      memset(seq_cmp,'\0',MAXS);
      nid=0;
      for (j=0;j<template[kk][k].main.length;j++)
      {
        itmp=template[kk][k].main.element[j]-1;
        for (jj=0;jj<MAXE;jj++){count[jj]=0;}
        for (jj=0;jj<template[kk][k].size;jj++){
            seq_tmp[jj][j]=traduce[template[kk][k].sequence[jj].element[itmp]];
            itmpir[jj][j]=template[kk][k].sequence[jj].position[itmp];
            count[template[kk][k].sequence[jj].element[itmp]]++;
        }
        jtmp=1000;
        for (jj=0;jj<MAXE;jj++){
             count[jj]=100*count[jj]/template[kk][k].size;
             if (ConsensusID(count[jj],template[kk][k].size)==1){jtmp=jj;break;}
             //if (count[jj]>50){jtmp=jj;break;}
        }
        if (jtmp < 1000){
         if (ConsensusID(count[jtmp],template[kk][k].size)==0){strncpy(&seq_cmp[j],".",1);}else{seq_cmp[j]=traduce[jtmp];nid++;}
        }else{
         strncpy(&seq_cmp[j],".",1);
        }
      }
      length=0;
      for (jj=0;jj<template[kk][k].size;jj++){memset(word_tmp[jj],'\0',MAXLINE);}
      memset(word_cmp,'\0',MAXLINE);
      for (j=0;j<template[kk][k].main.length;j++)
      {
            if (length==MAXLINE){
               fprintf(OUT,"\n");
               printf("\n");
               for (jj=0;jj<template[kk][k].size;jj++){
                 fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTmp[jj],itmpir[jj][j-length],word_tmp[jj],itmpir[jj][j-1]);
                 printf("%s %5d *->   %s <-*%5d \n",NameTmp[jj],itmpir[jj][j-length],word_tmp[jj],itmpir[jj][j-1]);
               }
               fprintf(OUT,"CONSENSUS_____             %s         \n",word_cmp);
               printf("CONSENSUS_____             %s         \n",word_cmp);
               length=0;
               for (jj=0;jj<template[kk][k].size;jj++){memset(word_tmp[jj],'\0',MAXLINE);}
               memset(word_cmp,'\0',MAXLINE);
            }      
            for (jj=0;jj<template[kk][k].size;jj++){ strncat(word_tmp[jj],&seq_tmp[jj][j],1); }
            strncat(word_cmp,&seq_cmp[j],1);
            length++;
      } 
      fprintf(OUT,"\n");
      printf("\n");
      for (jj=0;jj<template[kk][k].size;jj++){ 
        fprintf(OUT,"%s %5d *->   %s <-*%5d \n",NameTmp[jj],itmpir[jj][j-length],word_tmp[jj],itmpir[jj][j-1]);
        printf("%s %5d *->   %s <-*%5d \n",NameTmp[jj],itmpir[jj][j-length],word_tmp[jj],itmpir[jj][j-1]);
      }
      fprintf(OUT,"CONSENSUS_____             %s         \n",word_cmp);
      printf("CONSENSUS_____             %s         \n",word_cmp);
      
 

    free_cmatrix(NameTmp,0,template[kk][k].size,0,MAXTI);
    free_cvector(NameTmp0,0,MAXTI);

    free_cmatrix(seq_tmp,0,template[kk][k].size,0,MAXS);
    free_cmatrix(word_tmp,0,template[kk][k].size,0,MAXLINE);
    free_cvector(word_cmp,0,template[kk][k].size,0,MAXLINE);
    free_imatrix(itmpir,0,template[kk][k].size,0,MAXS);
    free_cvector(seq_cmp,0,MAXS);


}

