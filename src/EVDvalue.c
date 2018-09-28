#include  "inc.h"

float EVDevalue(parameter,score,size)
evd    parameter;
float  score,size;
{
  float  evalue;
  if (parameter.lambda * score > MAXEVD) {evalue=0.0;}
  else if (-1.0 * parameter.lambda * score > MAXEVD) { evalue= -1.0 * parameter.lambda * score;}
  else { evalue = parameter.K * size * exp( -1.0 * parameter.lambda * score);}
  return evalue;
}

float EVDpvalue(parameter,score,size)
evd    parameter;
float  score,size;
{
  float  pvalue,evalue;
  float  EVDevalue();
  evalue=EVDevalue(parameter,score,size);
  if (evalue>MAXEVD)    {pvalue= 1.0;}
  else if (evalue<1e-7) {pvalue=evalue;}
  else                  {pvalue= 1.0 - exp (-1.0 * evalue);}
  return pvalue;
}

