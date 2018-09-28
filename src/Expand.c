#include  "inc.h"

profile  Expand(p,a)
profile  p;
sequence a;
{
 profile r;
 profile InitProfile();
 int     i,j;

 r = InitProfile(a.title,a.length);
 r.size=p.size;
 r.main=a;
 r.cluster_dimension=p.cluster_dimension;
 for (i=0;i<p.cluster_dimension;i++) {r.cluster[i]=p.cluster[i];}
 for (i=0;i<p.size;i++) {
  strcpy(r.sequence[i].title,p.sequence[i].title);
  r.sequence[i].length=a.length;
  for (j=0;j<a.length;j++){
   if (a.element[j]>0){r.sequence[i].element[j]=p.sequence[i].element[a.element[j]-1];
                       r.sequence[i].position[j]=p.sequence[i].position[a.element[j]-1];}
   else               {r.sequence[i].element[j]=0;r.sequence[i].position[j]=r.sequence[i].position[j-1];}
  }
 }
 for (j=0;j<a.length;j++){
  if (a.element[j]>0){ r.profile[j]=p.profile[a.element[j]-1];}
  else               { r.profile[j].feature[0] ='-';
                       r.profile[j].feature[1] ='-';
                       r.profile[j].feature[2] ='-';
                       r.profile[j].feature[3] ='-';
                       r.profile[j].feature[4] ='-';
                     }
 }

 return  r;
}
