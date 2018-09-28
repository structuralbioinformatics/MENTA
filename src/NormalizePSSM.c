#include  "inc.h"

void NormalizePSSM(lena,lenb,m,n)
int     lena,lenb;
float **m,**n;
{
  int   ii,jj;
  float mu,mu2,m2u,sigma,minimum,maximum;
   

  minimum=MAXL;
  maximum=0.0;

  m2u=0.0;
  mu =0.0;
  for (ii=1;ii<=lena;ii++){
  for (jj=1;jj<=lenb;jj++){
     mu  += m[jj][ii];
     m2u += m[jj][ii] * m[jj][ii];
  }}

  mu   = mu / (lena*lenb) ;
  mu2  = mu * mu;
  m2u  = m2u / (lena*lenb) ;
  sigma= sqrtf ( fabs(m2u - mu2) );
  if (sigma < 1.0e-5) {sigma= 1.0;}

  for (ii=1;ii<=lena;ii++){ for (jj=1;jj<=lenb;jj++){ n[jj][ii] = 0.0; }}

  for (ii=1;ii<=lena;ii++){
  for (jj=1;jj<=lenb;jj++){
     n[jj][ii] = (m[jj][ii] - mu) / sigma;
     if (n[jj][ii]>0) n[jj][ii]  =  n[jj][ii] + 1.0;
     if (n[jj][ii]<0) n[jj][ii]  =  n[jj][ii] - 1.0;
  }}
  for (ii=0;ii<=lena;ii++){ n[0][ii]=0.0; }
  for (ii=0;ii<=lenb;ii++){ n[ii][0]=0.0; }

  for (ii=0;ii<=lena;ii++){
  for (jj=0;jj<=lenb;jj++){
    if ( n[jj][ii]>0 && minimum>n[jj][ii]) minimum=n[jj][ii];
    if ( n[jj][ii]>maximum)maximum=n[jj][ii];
  }}

  printf("\tAverage:\t%e\n",mu);
  printf("\tAverage²:\t%e\n",m2u);
  printf("\tSigma:  \t%e\n",sigma);
  printf("\tMaximum:\t%e\n",maximum);
  printf("\tMinimum:\t%e\n",minimum);

}

