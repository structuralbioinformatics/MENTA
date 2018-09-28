#include  "inc.h"

void  Align4Profile(a,template,target,result,matrix,traduce,n,number)
alignment a;
profile   *template;
profile   *target;
profile   *result;
float    **matrix;
char       traduce[MAXE];
int        n,number;
{
 int        i,k,j;
 float      ident,homol,simil[2];
 profile    r,p,q,dummy_tmp,dummy_tgt;
 profile    Expand();
 void       PushSeq(),Align2Similarity();
 profile    InitProfile(); 

 for (k=0;k<n;k++){
  j=0;
  dummy_tmp=template[k];
  dummy_tgt=target[k];
  p = Expand(dummy_tmp,a.template);
  q = Expand(dummy_tgt,a.target);
  r = InitProfile(a.title,a.length,number);
  r.score=a.score;
  r.evalue=a.evalue;
  r.pvalue=a.pvalue;
  Align2Similarity(a,dummy_tmp,dummy_tgt,matrix,traduce,simil);
  ident=simil[0];
  homol=simil[1];
  r.ident=ident;
  r.homol=homol;
  for (i=0;i<MAXS;i++){ r.profile[i]=a.profile[i];}
  for (i=0;i<p.size;i++) PushSeq( p.sequence[i],&r);
  for (i=0;i<p.cluster_dimension;i++) {r.cluster[j]=p.cluster[i];j++; r.cluster_dimension=j;}
  for (i=0;i<q.size;i++) PushSeq( q.sequence[i],&r);
  for (i=0;i<q.cluster_dimension;i++) {r.cluster[j]=q.cluster[i];j++; r.cluster_dimension=j;}
  result[k]=r;
 }

}
