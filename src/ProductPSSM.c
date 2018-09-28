#include  "inc.h"

void  ProductPSSM(lena,lenb,lenc,a,b,product)
int   lena,lenb,lenc;
float **a,**b,**product;
{

  int  i,j,k;

  for (j=1;j<=lenb;j++){
  for (k=1;k<=lenc;k++){
    product[j][k]=0.0;
    for (i=1;i<=lena;i++){ product[j][k] += a[j][i] * b[i][k];
  }}

  for (ii=0;ii<=lenb;ii++){ product[0][ii]=0.0; }
  for (ii=0;ii<=lenc;ii++){ product[ii][0]=0.0; }


}
