#include "inc.h"

float  CreateMatrix4PAIR(pssm,matrix,g,ge,p,pe,lms,method,psic_flag,tolerance,setpro,setmat,checkpro,seq,template,target,weight,weight_method,data,cte)
float  **pssm,**matrix,*g,*ge,*p,*pe;
float   *tolerance;
int      lms,method,setpro,setmat,checkpro,psic_flag;
char     seq[MAXP][MAXQ][MAXS];
profile *template,*target;
float    weight[MAXM],weight_method[MAXMETH],data[2][MAXM][MAXE][MAXE],cte[MAXP];
{

 int        mth,n,jk,ik,ii,jj,k,lena,lenb;
 float      maxim,maximum,ctpm,tolerance_pair[MAXTOL],**pssm_method,**pssm_normal,cte_pair[MAXP],weight_pair[MAXM];
 psic     **PSICtemplate,**PSICtarget;
 float    **qtemplate,**qtarget;
 alignment *result;
 sequence   Stemplate,Starget;
 float     *pb,*qb,*lambda1,**score;
 evd        evd_param;

 psic     **MEntPSIC();
 alignment *avector();
 float      CreateMatrix4METHOD(),PSSM2Matrix(),lambda(); 
 int        LocalDyna();
 float     *fvector(),**Fmatrix(),ffmax();
 void       PrintPSSM(),free_Fmatrix(),free_PSICmatrix(),free_avector(),free_fvector();

   lena=template[0].main.length;
   lenb=target[0].main.length;
   Stemplate=template[0].main;
   Starget  =target[0].main;

   qtemplate   =Fmatrix(0,MAXP+MAXM,0,MAXE);
   qtarget     =Fmatrix(0,MAXP+MAXM,0,MAXE);
   PSICtemplate=MEntPSIC(qtemplate,setpro,setmat,psic_flag,template,data,cte,weight);
   PSICtarget  =MEntPSIC(qtarget,setpro,setmat,psic_flag,target,data,cte,weight);

   evd_param.lambda=1.0;
   evd_param.K     =1.0;
   evd_param.sigma =1.0; 
   evd_param.Rc    =1;
   evd_param.Sc    =1.0;
   evd_param.min   =0.0;

   maxim=0.0;
   ctpm=0.0;
   for  (jk=0;jk<MAXMETH;jk++){ ctpm +=weight_method[jk];}
   for  (jk=0;jk<MAXMETH;jk++){ weight_method[jk] = weight_method[jk]/ctpm;}
   for  (jk=0;jk<=lenb; jk++){ for (ik=0;ik<=lena; ik++){ pssm[jk][ik]=0.0;}}


   printf("\n**************************\n");
   printf("Evaluate LAMBDA for bit-scores\n\n");

   lambda1=fvector(0,setpro+setmat);
   score=Fmatrix(0,MAXE,0,MAXE);
   pb=fvector(0,MAXE);
   qb=fvector(0,MAXE);
   
   lambda1[0]=1.0;
   for (n=1;n<setpro;n++){
   if (cte[n-1] != 0) {
    pb[0]=qb[0]=pb[1]=qb[1]=0.0;
    for (ii=2;ii<MAXE;ii++){
       score[0][ii]=score[ii][0]=score[1][ii]=score[ii][1]=0.0;
       pb[ii]=qtemplate[n][ii];
       qb[ii]=qtarget[n][ii];
       for (jj=2;jj<MAXE;jj++){score[ii][jj]=data[1][n-1][ii][jj];}
    }
    lambda1[n]=lambda(pb,qb,score,MAXE,MAXE);
    printf("Lambda Data[%d] = %e\n",n,lambda1[n]);
    if (lambda1[n] <=0.0) {cte[n-1]=0.0;lambda1[n]=1.0;}
   }}
   for (n=setpro;n<setpro+setmat;n++){
   k=n-setpro;
   if (weight[k]!=0){      
    pb[0]=qb[0]=pb[1]=qb[1]=0.0;
    for (ii=2;ii<MAXE;ii++){
       score[0][ii]=score[ii][0]=score[1][ii]=score[ii][1]=0.0;
       pb[ii]=qtemplate[0][ii];
       qb[ii]=qtarget[0][ii];
       for (jj=2;jj<MAXE;jj++){score[ii][jj]=data[0][k][ii][jj];}
    }
    lambda1[n]=lambda(pb,qb,score,MAXE,MAXE);
    printf("Lambda Data[%d] = %e\n",n,lambda1[n]);
    if (lambda1[n] <=0.0) {weight[k]=0.0;lambda1[n]=1.0;}
   }}   
 



   if (method==0){
     for (mth=1;mth<MAXMETH;mth++){ 
      if (weight_method[mth-1]>0){
       pssm_method  =Fmatrix(0,lenb+1,0,lena+1);
       result=avector(0,DIM);
       for  (jk=0;jk<=lenb; jk++){ for (ik=0;ik<=lena; ik++){ pssm_method[jk][ik]=0.0;}}
       printf("\n**************************\n"); 
       printf("\nMETHOD: %5d \n",mth);
       for(ik=0;ik<MAXTOL; ik++)tolerance_pair[ik]=tolerance[ik];
       for(ik=0;ik<MAXP;ik++)cte_pair[ik]=cte[ik];
       for(ik=0;ik<MAXM;ik++)weight_pair[ik]=weight[ik];
       maxim=0.0;
       maxim=CreateMatrix4METHOD(matrix,g,ge,p,pe,mth,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight_pair,data,cte_pair,lambda1,tolerance_pair); 
       if (lms >= 2){ maximum= (maxim)<=MAXV?(maxim):MAXV;}
       else         { maximum= (maxim/2)<=MAXV?(maxim/2):MAXV;}
       printf("Highest Score[%e]: %e \n",maxim,maximum); 
       //tolerance_pair[0]=tolerance[0];
       //tolerance_pair[1]=tolerance[1];
       //tolerance_pair[8]=tolerance[8];
       //tolerance_pair[9]=tolerance[9];
       for (jk=0;jk<MAXTOL;jk++){tolerance_pair[jk]=tolerance[jk];}
       if (tolerance[7]<= 0.0)   {tolerance_pair[6]= 500.0 * maximum;}else{ tolerance_pair[6]=tolerance[7];}
       if (tolerance[5]<= 0.0)   {tolerance_pair[2]= 500.0;}else{ tolerance_pair[2]=tolerance[5];}
       if (tolerance[3]>-1.0e+6) {tolerance_pair[3]=ffmax(0.0,tolerance[3]);}else{tolerance_pair[3]=tolerance[3];}
       tolerance_pair[4]=1.0; 
       if (mth == 5 || mth == 6 || mth == 7 || mth == 8 || mth == 13) tolerance_pair[4]=3.0;
       n=LocalDyna(matrix,Stemplate,Starget,g,ge,p,pe,lms,evd_param,maximum,tolerance_pair,pssm_method,result);
       printf("Local alignments: %d\n",n);
       if (n>0){printf("Weight of Scoring: %e\n",weight_method[mth-1]/n);}else{printf("Weight of Scoring: NULL\n");}
       if (n>0)for(jk=0;jk<=lenb; jk++)for(ik=0;ik<=lena; ik++) pssm[jk][ik]=pssm[jk][ik]+(weight_method[mth-1]/n)*pssm_method[jk][ik];
       free_avector(result,0,DIM);
       free_Fmatrix(pssm_method,0,lenb+1,0,lena+1);
      }
     }
     pssm_normal=Fmatrix(0,lenb+1,0,lena+1);
/*     PrintPSSM(lena,lenb,pssm);*/
     printf("\n**************************\n"); 
     printf("Normalize Alignement Matrices by Z-score\n");
     NormalizePSSM(lena,lenb,pssm,pssm_normal);
/*     PrintPSSM(lena,lenb,pssm_normal);*/
     maxim=0.0;
     maxim=PSSM2Matrix(pssm_normal,matrix,g,ge,p,pe,setpro,setmat,checkpro,template,target,qtemplate,PSICtemplate,PSICtarget,qtarget,weight,data,cte,lambda1,tolerance);
     free_Fmatrix(pssm_normal,0,lenb+1,0,lena+1);
   }else{
      printf("\n**************************\n"); 
      printf("\nMETHOD: %5d\n",method);
      maxim=0.0;
      (maxim)=CreateMatrix4METHOD(matrix,g,ge,p,pe,method,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance); 
   }


      free_Fmatrix(qtemplate,0,MAXP+MAXM,0,MAXE);
      free_Fmatrix(qtarget,0,MAXP+MAXM,0,MAXE);
      free_PSICmatrix(PSICtemplate,0,setpro+setmat,0,lena);
      free_PSICmatrix(PSICtarget,0,setpro+setmat,0,lenb);
      free_fvector(lambda1,0,setpro+setmat);
      free_Fmatrix(score,0,MAXE,0,MAXE);
      free_fvector(pb,0,MAXE);
      free_fvector(qb,0,MAXE);


      return  maxim ;

}
