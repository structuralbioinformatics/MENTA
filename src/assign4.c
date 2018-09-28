#include  "inc.h"

void  Assign4(seq,title,ii,jj,traduce,setpro,n_set,id_set,set,template,target)
char      title[MAXP][MAXQ][50];
char      seq[MAXP][MAXQ][MAXS];
int       ii,jj,setpro;
char      traduce[MAXE];
int       *n_set,*id_set,**set;
profile  *template,*target;
{
 int      i,j,ji,k,ipos,lena,lenb;
 char     title_profile[100+MAXS],set_of_sequences[MAXS];
 sequence *a,*b;
 sequence *svector();
 void     free_svector();
 profile  Seq2Profile();
 
 lena=strlen(seq[0][set[ii][0]]);
 lenb=strlen(seq[0][set[jj][0]]);
 a=svector(0,n_set[ii]);
 b=svector(0,n_set[jj]);


 for (k=0;k<setpro;k++){
  for (ji=0;ji<n_set[ii];ji++){
    strcpy(a[ji].title,title[k][set[ii][ji]]);
    ipos=0;
    a[ji].length=lena;
    for (j=0;j<lena;j++){
      a[ji].element[j]=100;
      for (i=0;i<MAXE;i++){
        if (traduce[i]==seq[k][set[ii][ji]][j]){a[ji].element[j]=i; break;}
      }
      if (a[ji].element[j]==100) {printf("Target Data: %d [%d][%d] Property: %d Position: %d Seq: %s Title: %s\n",set[ii][ji],ii,ji,k,j,seq[k][set[ii][ji]][j],a[ji].title);
                                  nrerror("Undefined residue code in template sequence (Check length)");}
      if (strncmp(&traduce[a[ji].element[j]],"-",1)){ ipos++;a[ji].position[j]=ipos;}else{a[ji].position[j]=ipos;}
    }
  }
  for (ji=0;ji<n_set[jj];ji++){
    strcpy(b[ji].title,title[k][set[jj][ji]]);
    ipos=0;
    b[ji].length=lenb;
    for (j=0;j<lenb;j++){
      b[ji].element[j]=100;
      for (i=0;i<MAXE;i++){
        if (traduce[i]==seq[k][set[jj][ji]][j]){b[ji].element[j]=i; break;}
      }
      if (b[ji].element[j]==100) {printf("Template Data: %d  [%d][%d] Property: %d Position: %d Seq: %s Title: %s\n",set[jj][ji],jj,ji,k,j,seq[k][set[jj][ji]][j],b[ji].title);
                                  nrerror("Undefined residue code in target sequence (Check length)");}
      if (strncmp(&traduce[b[ji].element[j]],"-",1)){ipos++;b[ji].position[j]=ipos;}else{b[ji].position[j]=ipos;}
    }
  }

  memset(set_of_sequences,'\0',MAXS);
  for (ji=0;ji<n_set[ii];ji++){if ((strlen(set_of_sequences) + 10)<MAXS){ sprintf(set_of_sequences,"%s%d-",set_of_sequences,set[ii][ji]);}}
  sprintf(title_profile,"Profile %d containing %d sequences: %s \n",id_set[ii],n_set[ii],set_of_sequences);
  template[k]= Seq2Profile(title_profile,a,n_set[ii]);
  memset(set_of_sequences,'\0',MAXS);
  for (ji=0;ji<n_set[jj];ji++){if ((strlen(set_of_sequences) + 10)<MAXS){sprintf(set_of_sequences,"%s%d-",set_of_sequences,set[jj][ji]);}}
  sprintf(title_profile,"Profile %d containing %d sequences: %s \n",id_set[jj],n_set[jj],set_of_sequences);
  target[k]  = Seq2Profile(title_profile,b,n_set[jj]);

 }
 free_svector(a,0,n_set[ii]);
 free_svector(b,0,n_set[jj]);
}
