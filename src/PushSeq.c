#include  "inc.h"

void  PushSeq(a,p)
sequence a;
profile  *p;
{
 int   length,size,i,j,skip;
 void  nrerror();
 char  error_size[MAXS];
 char  error_length[MAXS];

 length = p[0].main.length;
 size   = p[0].size;
 sprintf (error_size,"Too Many Sequences on a Profile: %s",p[0].main.title);
 sprintf (error_length,"Incompatible Sequence length for the Profile: %s",p[0].main.title);
 if ( size   >=  MAXFS     ) nrerror(error_size);
 if ( length !=  a.length ) nrerror(error_length);
 p[0].sequence[size]=a;
 p[0].size=size+1;
 p[0].align=0;
 for (j=0;j<length;j++){
   skip=0;
   for ( i=0;i<p[0].size;i++) {if (p[0].sequence[i].element[j] <=1 ) {skip=1;break;}}
   if (skip==0){p[0].align++;}
 }
 return;

}
