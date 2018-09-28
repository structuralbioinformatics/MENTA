#include  "inc.h"

 float lambda(p,q,s,n,m)
 float *p,*q;
 float **s;
 int   n,m;
 {
   float a,b,fa,fb,fb0,fa0,sign,sum,suma,sumb,diff;
   int   i,j,k,kk,gon;
   float res,x;

  if (EVAL == 3) return 1.0;

  kk=0;
  gon=0;
  for (i=0;i<n;i++) for (j=0;j<m;j++) if (p[i]>0.0 && q[j]>0.0 && s[i][j] >0) kk++;
  if (kk==0){gon=1;}
  kk=0;
  if (gon==1){
   for (i=0;i<n;i++) for (j=0;j<m;j++) if (s[i][j]>0){ kk++;if (p[i]==0.0)p[i]=1.0e-8;if (q[j]==0.0)q[j]=1.0e-8;}
   if (kk==0){res=0.0;printf("Escape Lambda (no positive scores)\n");return res;}
  }

  kk=0;
  gon=0;
  for (i=0;i<n;i++) for (j=0;j<m;j++) if (p[i]>0.0 && q[j]>0.0 && s[i][j] <0) kk++;
  if (kk==0){gon=1;}
  kk=0;
  if (gon==1){
   for (i=0;i<n;i++) for (j=0;j<m;j++) if (s[i][j]<0){ kk++;if (p[i]==0.0)p[i]=1.0e-6;if (q[j]==0.0)q[j]=1.0e-6;}
   if (kk==0){res=0.0;printf("Escape Lambda (no negative scores)\n");return res;}
  }

  a=0.5;
  b=1.0;
  sum=0.0;
  for (i=0;i<n;i++){ for (j=0;j<m;j++){
             if (p[i]>0.0 && q[j]>0.0){sum += p[i]*q[j]*exp(b * s[i][j]);  }}}
  fb=1.0-sum;
  fb0=fb;
  k=0;
  kk=0;
  sign=1.0;
  while(fb>0.0 && k<1.0e+2){
   b+=1.0;
   k++;
   kk++;
   sum=0.0;
   for (i=0;i<n;i++){ for (j=0;j<m;j++){
             if (p[i]>0.0 && q[j]>0.0){sum += p[i]*q[j]*exp(b * s[i][j]); }}}
   if (kk==10 ){kk=0;if (fb0<fb) sign=-1.0;fb0=fb;printf("Lambda FB(%e)= %e\n",b,fb);}
   fb=1.0-sum;
  }
  if (fb>0.0){printf("No convergence for Lambda FB(%e)= %e\n",b,fb);}

  sum=0.0;
  for (i=0;i<n;i++){ for (j=0;j<m;j++){
             if (p[i]>0.0 && q[j]>0.0){sum += p[i]*q[j]*exp(a * s[i][j]); }}}
  fa=1.0-sum;
  fa0=fa;
  k=0;
  kk=0;
  sign=1.0;
  while(fa<0.0 && k<1.0e+2){
   a/=1.5;
   k++;
   kk++;
   sum=0.0;
   for (i=0;i<n;i++){ for (j=0;j<m;j++){
             if (p[i]>0.0 && q[j]>0.0){sum += p[i]*q[j]*exp(a * s[i][j]); }}}
   if (kk==10 ){kk=0;if(fa0>fa) sign=-1.0*sign;fa0=fa;printf("Lambda FA(%e)= %e\n",a,fa);}
   if ( a <1.0e-12)break;
   fa=1.0-sum;
  }

  if (fa<0.0){printf("No convergence for Lambda FA(%f)= %f\n",a,fa);res=0.0;return res;}

  if (a>b) {a=b;b=a;}
  k=0;
  res=a+(b-a)/2.0;
  sum=0.0;
  for (i=0;i<n;i++){ for (j=0;j<m;j++){
             if (p[i]>0.0 && q[j]>0.0){sum += p[i]*q[j]*exp(res * s[i][j]); }}}
  diff=1.0-sum;


  while (fabs(diff) >= 1.0e-8 && k < 100){
    x=fabs(b-a)/2.0;
    if (fa*fb<0){
     if (diff*fb > 0.0 ){
       b=res;
       res=a+x;
     }
     if (diff*fa>0.0){
       a=res;
       res=b-x;
     }
     if (fabs(b-a)<1.0e-8){res=a;k=100;}
    }else{
      if (fb>0.0){
            kk=0;
            while(fb>0.0 && kk<100){
               b+=1.0;
               kk++;
               sum=0.0;
               for (i=0;i<n;i++){ for (j=0;j<m;j++){
                 if (p[i]>0.0 && q[j]>0.0){sum += p[i]*q[j]*exp(b * s[i][j]); }}}
               fb=1.0-sum;
            }
      }else if (fa<0){
            kk=0;
            while(fa<0.0 && kk<100){
             a+=0.5;
             kk++;
             sum=0.0;
             for (i=0;i<n;i++){ for (j=0;j<m;j++){
                       if (p[i]>0.0 && q[j]>0.0){sum += p[i]*q[j]*exp(a * s[i][j]); }}}
             fa=1.0-sum;
            }
      } 
      if (a>b) {a=b;b=a;}
      res=a+(b-a)/2.0;   
    }

    sumb=0.0;
    for (i=0;i<n;i++){ for (j=0;j<m;j++){
             if (p[i]>0.0 && q[j]>0.0){sumb += p[i]*q[j]*exp(b * s[i][j]); }}}
    fb=1.0-sumb;
    suma=0.0;
    for (i=0;i<n;i++){ for (j=0;j<m;j++){
             if (p[i]>0.0 && q[j]>0.0){suma += p[i]*q[j]*exp(a * s[i][j]); }}}
    fa=1.0-suma;
    sum=0.0;
    for (i=0;i<n;i++){ for (j=0;j<m;j++){
             if (p[i]>0.0 && q[j]>0.0){ sum += p[i]*q[j]*exp(res*s[i][j]); }}}
    diff=1.0-sum;
    k++;
  }
  if (res<=0.0)res=0.0;
  return res;

 }
