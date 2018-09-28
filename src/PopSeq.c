#include  "inc.h"

sequence  PopSeq(p,j)
profile  *p;
int      j;
{
 int       i,length;
 char      title[MAXTITLE];
 sequence  a;
 profile   *q;
 profile   *pvector();
 profile   InitProfile();
 void      free_pvector();
 
 length=p[0].main.length;
 strcpy(title,p[0].main.title);
 a = p[0].sequence[j];
 q = pvector(0,1);
 q[0] = InitProfile(title,length);
 for (i=0;i<p[0].size;i++)  if (i!=j) PushSeq(p[0].sequence[i],q);
 p[0] = q[0];
 free_pvector(q,0,1);

 return a;
}
