#include  "inc.h"

void SumPSSM(lena,lenb,a,b,sum)
int     lena,lenb;
float **a,**b,**sum;
{
  int   ii,jj;

  for (ii=1;ii<=lena;ii++){for (jj=1;jj<=lenb;jj++){sum[jj][ii] =0.0;}}
  for (ii=1;ii<=lena;ii++){
  for (jj=1;jj<=lenb;jj++){
    sum[jj][ii] = a[jj][ii]+  b[jj][ii];
  }}
  for (ii=0;ii<=lena;ii++){ sum[0][ii]=0.0; }
  for (ii=0;ii<=lenb;ii++){ sum[ii][0]=0.0; }

}

