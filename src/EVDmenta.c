#include   "inc.h"

evd EVDmenta(m,a,b,g,ge,p,pe,lms,limit,tolerance)
float         **m;
sequence      a;
sequence      b;
float          *g,*ge,*p,*pe;
int             lms;
float          limit;
float         *tolerance;
{
 float          lambda_evd,minimum;
 sequence      template,target;
 int             i,j,n,lena,lenb,sumR;
 evd             parameter,evd_island;
 int             ik,jk;
 float          *pb,*qb,**score;


 int             EVDisland();
 float       **Fmatrix();
 float         *fvector();
 void            free_Fmatrix();
 void            free_vector();
 float          lambda(),ffmax();



  if   (EVAL >= 2 ){

    parameter.lambda=1.0;
    parameter.K       =1.0;
    parameter.sigma   =1.0; 
    parameter.Rc      =1;
    parameter.Sc      =1.0;
    parameter.min     =0.0;

  }else{

    lena=a.length;
    lenb=b.length;
    pb=fvector(0,lena);
    qb=fvector(0,lenb);
    score=Fmatrix(0,lena,0,lenb);

    printf("Extreme Value Distribution Parameters\n");
      for (ik=0;ik<lena;ik++){
       pb[ik]=1.0/lena;
       for (jk=0;jk<lenb;jk++){
            qb[jk]=1.0/lenb;
            score[ik][jk] =   m[jk][ik] ;
         }
      }
    lambda_evd=lambda(pb,qb,score,lena,lenb);

    sumR=EVDisland(m,a,b,g,ge,p,pe,lms,&evd_island,limit,tolerance);

    minimum=ffmax(tolerance[3],evd_island.min);

    switch (EVAL){
      case 0:
               parameter.lambda=lambda_evd;
               parameter.K= sumR * exp ( lambda_evd * minimum) / (a.length * b.length);
               parameter.sigma=1.0;
               printf("EVD Pairs   (MENTA)= %d\n",sumR);
               printf("EVD Lambda (MENTA)= %e\n",parameter.lambda);
               printf("EVD K         (MENTA)= %e\n",parameter.K);
               break;
      case 1:
               parameter=evd_island;
               printf("EVD Lambda (ISLAND)= %e\n",parameter.lambda);
               printf("EVD K       (ISLAND)= %e\n",parameter.K);
               printf("EVD Sigma   (ISLAND)= %e\n",parameter.sigma);
               break;
      default:
               parameter.lambda=lambda_evd;
               parameter.K= sumR * exp ( lambda_evd * minimum ) / (a.length * b.length);
               parameter.sigma=1.0;
               printf("EVD Pairs   (MENTA)= %d\n",sumR);
               printf("EVD Lambda (MENTA)= %e\n",parameter.lambda);
               printf("EVD K         (MENTA)= %e\n",parameter.K);
               break;
      }

    free_Fmatrix(score,0,lena,0,lenb);
    free_fvector(pb,0,lena);
    free_fvector(qb,0,lenb);

  }

  return parameter;
}
