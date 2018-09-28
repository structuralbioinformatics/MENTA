#include  "inc.h"

float  extension(g,ipos,jext,limit)
float       *g,*limit;
int          ipos;
int          jext;
{
  float      ge;
  int        i;
  ge=0.0;
  for  (i=ipos-jext-1;i<ipos-1;i++){ge += g[i];}
  if  ( (ipos-jext-1) < 2 ) {ge=0.0;} 
  if  ( ipos < 1 ) {ge=0.0;} 
  if ( ge > limit[6]) {ge=0;}
  if ( ipos-jext-1 > limit[2]) {ge=0;}
  return     ge;
}
