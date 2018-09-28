#include  "inc.h"

void PushPSSM(n,lena,lenb,m,ma)
int     n,lena,lenb;
float **m;
float ***ma;
{
  int   i,j,k;

  for (k=0;k<n;k++){
  for (i=0;i<lena;i++){
  for (j=0;j<lenb;j++){
     ma[n][j][i]=0.0;
  }}}

}
