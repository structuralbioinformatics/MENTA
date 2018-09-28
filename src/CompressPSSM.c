#include  "inc.h"

void CompressPSSM(lena,lenb,m,n)
int     lena,lenb;
float **m,**n;
{
  int   ii,jj,k;
  float max,ave;

  k=0;
  max=1.0;
  ave=0.0;

  for (ii=1;ii<=lena;ii++){
  for (jj=1;jj<=lenb;jj++){
   if (m[jj][ii]!=0){
      k++;
      max = max >= abs(m[jj][ii]) ? max:abs(m[jj][ii]) ;
      ave += abs(m[jj][ii]);
   }
  }}

  ave = ave/k;

  for (ii=1;ii<=lena;ii++){
  for (jj=1;jj<=lenb;jj++){
    if (ave > 1.0){ 
     n[jj][ii] = m[jj][ii] / ave ;
    }else{
     n[jj][ii] = m[jj][ii] / max ;
    }
  }}
  for (ii=0;ii<=lena;ii++){ n[0][ii]=0.0; }
  for (ii=0;ii<=lenb;ii++){ n[ii][0]=0.0; }

}

