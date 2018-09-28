#include "inc.h"

psic  **MEntPSIC(q,setpro,setmat,pattern,data)
int      setpro,setmat;
float  **q;
profile *pattern;
float   data[2][MAXM][MAXE][MAXE];
{
 psic  **Kroenecker,**PSICpattern;
 int     i,j,k,n,ii,jj,kk,len,size;
 float   Neff,Nobs,EqPos,product,ratio,z,g,alfa,beta,alfaAA,betaAA;
 float   frequency[MAXE]; 
     
 float   Bernoulli();
 psic  **PSICmatrix();
 void    free_PSICmatrix();


 len =pattern[0].main.length;
 size=pattern[0].size;
 Kroenecker  =PSICmatrix(0,size,0,len);
 PSICpattern =PSICmatrix(0,setpro+setmat,0,len);

 printf ("\nPseudo Counts (PSIC)\n");
 for (n=0;n<setpro;n++){
/*----------------------------------------------------------------------------------*/
  printf ("Calculating Pseudo Counts (PSIC) for property %d \n",n);
  for (k=0;k<MAXE;k++) {frequency[k]=0.0;}
  for (i=0;i<size;i++){
   for (j=0;j<len;j++){
       for (k=0;k<MAXE;k++) {Kroenecker[i][j].element[k] =0.0;}
       Kroenecker[i][j].element[pattern[n].sequence[i].element[j]]=1.0;
       frequency[pattern[n].sequence[i].element[j]] += 1.0;
   }
  }
  z=0.0;
  for (k=2;k<MAXE;k++){ z+=frequency[k];}
  if (size>=1 && len >=1 && z>0) {
     for (k=0;k<2;k++){frequency[k] = frequency[k] /(size*len);};
     for (k=2;k<MAXE;k++) {frequency[k] = frequency[k] / z;}
  }
  for (k=0;k<MAXE;k++) {q[n][k]=frequency[k];}
  alfa=0.0;
  beta=0.0;
  for (k=0;k<MAXE;k++){
   for (j=0;j<len;j++){
     Nobs =0.0;
     EqPos=0.0;
     for (i=0;i<size;i++){ Nobs+= Kroenecker[i][j].element[k]; }
     if (Nobs>0){alfa+=1.0;}
     for (jj=0;jj<len;jj++){
         if (jj!=j){ 
                    for (kk=0;kk<MAXE;kk++){
                        product=1.0;
                        for (ii=0;ii<size;ii++){ product *= Kroenecker[ii][j].element[k]*Kroenecker[ii][jj].element[kk];}
                        EqPos+=product;
                    }
         }
     }
     if (pattern[0].align > 0 && EqPos > 1) { 
       ratio=EqPos/pattern[0].align;
       Neff=Bernoulli(ratio,frequency,Nobs);
     }else{
       Neff=Nobs;
     }
     PSICpattern[n][j].element[k]=Neff;
   }
   if (Nobs>0 && k>1 ){beta+=1.0;}
  }
  alfa=  (alfa/((float) len)) - 1.0;
  beta=  beta/2.0;
  for (j=0;j<len;j++){
      z=0.3;
      for (k=0;k<MAXE;k++){ z+=PSICpattern[n][j].element[k];}
      for (k=0;k<MAXE;k++){ 
               if (PSICpattern[n][j].element[k]>0.0) {PSICpattern[n][j].frequency[k]= PSICpattern[n][j].element[k]/z;}
               else                                  {PSICpattern[n][j].frequency[k]= frequency[k]*0.3/z;}
      }
  }

  if (n>0){ 
   for (j=0;j<len;j++){
     for (k=2;k<MAXE;k++){
       g = 0.0;
       for (kk=2;kk<MAXE;kk++){
          g +=   PSICpattern[n][j].frequency[kk] * frequency[k] * ( ( log(data[1][n-1][k][kk])/log(2.0) ) - 1.0 ) ;
       }
       PSICpattern[n][j].targetQ[k]= (alfa*PSICpattern[n][j].frequency[k] + beta * g ) /(alfa + beta) ;
     }
     for (k=0;k<2;k++){PSICpattern[n][j].targetQ[k]=PSICpattern[n][j].frequency[k];}
   } 
  }else{
    for (j=0;j<len;j++){for (k=0;k<MAXE;k++){PSICpattern[n][j].targetQ[k]=0.0;}}
    alfaAA=alfa;
    betaAA=beta;
  }
/*----------------------------------------------------------------------------------*/
 }
 for (jj=setpro;jj<setmat+setpro;jj++){
  printf ("Calculating Pseudo Counts (PSIC) for PSSM %d \n",jj-setpro);
    for (j=0;j<len;j++){
       for (k=2;k<MAXE;k++){
         g = 0.0;
         for (kk=2;kk<MAXE;kk++){
          g +=   PSICpattern[0][j].frequency[kk] * q[0][k] * ( ( log(data[0][jj-setpro][k][kk])/log(2.0) ) - 1.0 ) ;
         }
         PSICpattern[jj][j].targetQ[k]= (alfaAA*PSICpattern[0][j].frequency[k] + betaAA * g ) /(alfaAA + betaAA) ;
         PSICpattern[jj][j].frequency[k]=PSICpattern[0][j].frequency[k];
         PSICpattern[jj][j].element[k]=PSICpattern[0][j].element[k];
       }
       for (k=0;k<2;k++){PSICpattern[jj][j].targetQ[k]=PSICpattern[0][j].frequency[k];
                         PSICpattern[jj][j].frequency[k]=PSICpattern[0][j].frequency[k];
                         PSICpattern[jj][j].element[k]=PSICpattern[0][j].element[k];
                        }
     }
 }
 
 free_PSICmatrix(Kroenecker,0,size,0,len);
 return PSICpattern; 

}

float Bernoulli(p,q,n)
float  p;
float  q[MAXE];
float  n;
{
  float result,sum,a,b,d;
  int   i,k;
  
  a=1.0;
  b=n;
  result=(b-a)/2.0;
  sum=0.0;
  for (i=0;i<MAXE;i++){ sum+=(float)pow((double)q[i],(double)result); }
  d=p-sum;
  k=0;
  while (abs(d) >= 1.0e-4 && k < 100){
    if (d<0.0){
       b=result;
       result=(b-a)/2.0;
    }else{
       a=result;
       result=(b-a)/2.0;
    }
    sum=0.0;
    for (i=0;i<MAXE;i++){ sum+=(float)pow((double)q[i],(double)result); }
    d=p-sum;
    k++;
  }

  return result;

}
