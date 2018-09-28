#include  "inc.h"

void  Assign(seq,title,ii,traduce,setpro,n_set,id_set,set,template,number)
char      title[MAXP][MAXQ][MAXTITLE];
char      seq[MAXP][MAXQ][MAXS];
int       ii,setpro;
char      traduce[MAXE];
int       *n_set,*id_set,**set;
profile  *template;
int       number;
{
 int      i,j,ji,k,ipos,lena,lenb;
 char     title_profile[MAXTITLE],set_of_sequences[MAXS];
 sequence *a;
 sequence *svector();
 void     free_svector();
 profile  Seq2Profile();
 
 lena=strlen(seq[0][set[ii][0]]);
 a=svector(0,n_set[ii]);


 for (k=0;k<setpro;k++){
  for (ji=0;ji<n_set[ii];ji++){
    //printf("Assign[%d][%d] Sequence Name %s\n",k,ji,title[0][set[ii][ji]]);
    strncpy(a[ji].title,title[k][set[ii][ji]],MAXTITLE);
    ipos=0;
    a[ji].length=lena;
    for (j=0;j<lena;j++){
      a[ji].element[j]=100;
      for (i=0;i<MAXE;i++){
        if (traduce[i]==seq[k][set[ii][ji]][j]){a[ji].element[j]=i; break;}
      }
      if (a[ji].element[j]==100) {printf("Sequence Data#: %d [%d][%d] Property: %d Position: %d AA: %c Title: %s\n",set[ii][ji],ii,ji,k,j,seq[k][set[ii][ji]][j],a[ji].title);
                                  nrerror("Undefined residue code in template sequence (Check length)");}
      if (strncmp(&traduce[a[ji].element[j]],"-",1)){ ipos++;a[ji].position[j]=ipos;}else{a[ji].position[j]=ipos;}
    }
  }

  memset(set_of_sequences,'\0',MAXS);
  for (ji=0;ji<n_set[ii];ji++){
     if (strlen(set_of_sequences )<MAXTITLE-200){ 
        sprintf(set_of_sequences,"%s%d-",set_of_sequences,set[ii][ji]);
     }else{
        nrerror("Too long title for Profile, increase MAXTITLE");
     }
   }
  sprintf(title_profile,"Profile_init %d containing %d sequences: %s \n",id_set[ii],n_set[ii],set_of_sequences);
  //printf("Assign Set[%d][%d] Profile[%d] containing %d sequences: %s \n",ii,k,id_set[ii],n_set[ii],set_of_sequences);
  template[k]= Seq2Profile(title_profile,a,n_set[ii],number);
  template[k].cluster[0]=ii;
  template[k].cluster_dimension=1;
  template[k].score=100.0;  
  template[k].evalue=0.0;  
  template[k].pvalue=0.0;  
  template[k].ident=100.0;
  template[k].homol=100.0;

 }

 free_svector(a,0,n_set[ii]);
}
