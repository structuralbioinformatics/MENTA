#include  "inc.h"

void SelectProfiles(dimension,number_of_profiles,setpro,tolerance,profiles)
int     dimension,number_of_profiles,setpro;
float   tolerance[MAXTOL];
profile **profiles;
{ 

 int   i,j,k,ii,jj,nn,dim;
 int   *grouped,*n,**group,skip[MAXPROF];
 void  gShell(),nrerror(),free_ivector(),free_imatrix();
 int   *ivector(),**imatrix();


   grouped=ivector(0,MAXPROF);
   n=ivector(0,MAXPROF);
   group=imatrix(0,MAXPROF,0,MAXPS);

   for (i=0;i<MAXPROF;i++){grouped[i]=-1;n[i]=0;skip[i]=0;}
   for (i=0;i<MAXPROF;i++){for (j=0;j<MAXPS;j++){group[i][j]=-1;}}

   
   for (ii=dimension;ii<number_of_profiles;ii++){
    dim=profiles[ii][0].cluster_dimension;
    if (grouped[ii]!=-1){
       n[ii]=n[grouped[ii]];
       for (jj=0;jj<n[ii];jj++){ group[ii][jj]=group[grouped[ii]][jj]; }
    }
    if (grouped[ii]==-1){
     grouped[ii]=ii;
     group[ii][n[ii]]=ii;
     n[ii]++;
     if (n[dim]>=MAXPS)nrerror("Too Many Sequences on SelectProfiles: Increase MAXPS");
     for (jj=ii+1;jj<number_of_profiles;jj++){
      dim=0;
      if (profiles[ii][0].cluster_dimension==profiles[jj][0].cluster_dimension){
       for (k=0;k<profiles[ii][0].cluster_dimension;k++){
       for (j=0;j<profiles[jj][0].cluster_dimension;j++){
          if (profiles[ii][0].cluster[k]==profiles[jj][0].cluster[j])dim++;
          }}
       if (dim==profiles[ii][0].cluster_dimension){
          group[grouped[ii]][n[grouped[ii]]]=jj;
          grouped[jj]=grouped[ii];
          n[grouped[ii]]++;
          }
      }
     }
    }
   }

   for (i=dimension;i<number_of_profiles;i++){
     if (profiles[i][0].action==1 && skip[i]==0){
        if (n[i]>tolerance[10]){ 
          gShell(i,n,group,profiles); 
          for (j=0;j<n[i];j++){skip[group[i][j]]=1;}
          for (j=(int)tolerance[10];j<n[i];j++){
            for (k=0;k<setpro;k++){profiles[group[i][j]][k].action=0;}
          }
        }
     }
   }
   
   free_ivector(grouped,0,MAXPROF);
   free_ivector(n,0,MAXPROF);
   free_imatrix(group,0,MAXPROF,0,MAXPS);
   
}
