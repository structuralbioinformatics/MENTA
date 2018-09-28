#include  "inc.h"

profile  Seq2Profile(title,a,d,number)
char      title[MAXTITLE];
sequence *a;
int       d,number;
{
 int       i,length;
 profile   p;
 profile   InitProfile();
 void      PushSeq();

 length = a[0].length;
 p = InitProfile(title,length,number);
 for (i=0;i<d;i++) if ( length == a[i].length ) PushSeq(a[i],&p);
 return  p;
}
