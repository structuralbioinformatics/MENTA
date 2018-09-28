#include  "inc.h"

void PushPSSM(n,lena,lenb,m,ma)
int     n,lena,lenb;
float **m;
float ***ma;
{
  int   i,j;

  for (i=0;i<=lena;i++){
  for (j=0;j<=lenb;j++){
     ma[n][j][i]=m[j][i];
  }}

}
