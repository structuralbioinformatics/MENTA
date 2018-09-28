#include  "inc.h"

void AccumPSSM(lena,lenb,a,b)
int     lena,lenb;
float **a,**b;
{
  int   ii,jj;

  for (ii=1;ii<=lena;ii++){
  for (jj=1;jj<=lenb;jj++){
    if (fabs(a[jj][ii]) > 1.0e-10) b[jj][ii] += a[jj][ii];
  }}
  for (ii=0;ii<=lena;ii++){ b[0][ii]=0.0; }
  for (ii=0;ii<=lenb;ii++){ b[ii][0]=0.0; }

}

