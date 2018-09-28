#include  "inc.h"

double  evalue(parameter,align,size)
evd          parameter;
alignment    align;
double       size;
{
 double  e_value,p_value;
 double  score;
 double  EVDevalue();
 double  pvalue();

 score=align.score;
 if (parameter.Rc < 10){
  p_value= pvalue(parameter,align,size);
  e_value= DATASIZE * p_value;
 }else{
  e_value=EVDevalue(parameter,score,size);
 }

 return e_value;

}

double  pvalue(parameter,align,size)
evd          parameter;
alignment    align;
double       size;
{
 double  p_value,factn,factk;
 double  score;
 double  EVDpvalue(),lfactorial();
 int     n,d,k;

 

 score=align.score;

 if (parameter.Rc < 10){
   n    =align.length;
   k    =(int) (n * align.ident);
   d    = n - k;
   factn=lfactorial(n,k);
   factk=lfactorial(k,1);
   p_value = exp (  factn -  factk + k * log (PMURZIN) + d * log ( 1.0 - PMURZIN ) );  
 }else{
   p_value=EVDpvalue(parameter,score,size);
 }

 return p_value;

}


double  EVDevalue(parameter,score,size)
evd    parameter;
double   score,size;
{
  double   evalue;
  if (parameter.lambda * score > MAXEVD) {evalue=0.0;}
  else if (-1.0 * parameter.lambda * score > MAXEVD) { evalue= -1.0 * parameter.lambda * score;}
  else { evalue = parameter.K * size * exp( -1.0 * parameter.lambda * score);}
  return evalue;
}

double  EVDpvalue(parameter,score,size)
evd      parameter;
double   score,size;
{
  double   pvalue,evalue;
  double   EVDevalue();
  evalue=EVDevalue(parameter,score,size);
  if (evalue>MAXEVD)    {pvalue= 1.0;}
  else if (evalue<1e-7) {pvalue=evalue;}
  else                  {pvalue= 1.0 - exp (-1.0 * evalue);}
  return pvalue;
}

double lfactorial(n,m)
int n,m;
{
  int    i;
  double z,fact;
  double lfactorial();


  if ( n<m || n<=0 || m<=0){z=0.0;}else{z=log((double)n);}

  if ( z > 0 ){
     i=n--;
     fact = z + lfactorial(i,m);
  }else{
     fact = z;
  }

  return fact;
}
