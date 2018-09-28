#include  "inc.h"

profile  InitProfile(a,l,n)
char    *a;
int      l,n;
{
 profile p;
 int     i;

 sprintf(p.main.title,"PROFILE[%d] ",n);
 //strncat(p.main.title,a,MAXTITLE-15);
 p.main.length=l;
 for (i=0;i<p.main.length;i++) p.main.element[i]=i+1;
 p.size=0;
 p.align=0;
 p.action=1;

 return p;
}
