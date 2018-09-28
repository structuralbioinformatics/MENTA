#include  "inc.h"

sequence  XtractSeq(p,j)
profile  p;
int      j;
{
 sequence a;
 a = p.sequence[j];
 return a;
}
