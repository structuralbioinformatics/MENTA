#include "inc.h"

float  CreateMatrix4METHOD(matrix,g,ge,p,pe,method,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight0,data,cte0,lambda1,tolerance)
float  **matrix,*g,*ge,*p,*pe,*lambda1,*tolerance;
int      method,setpro,setmat,checkpro;
char     seq[MAXP][MAXQ][MAXS];
profile *template,*target;
psic   **PSICtemplate,**PSICtarget;
float  **qtemplate,**qtarget;
float    weight0[MAXM],data[2][MAXM][MAXE][MAXE],cte0[MAXP];
{
 int lena,lenb,i,j,n,ii,jj,k;
 float maxim;
 float weight[MAXM],cte[MAXP];

 float CreateMatrix4raw();
 float CreateMatrix4SF();         
 float CreateMatrix4SQ();         
 float CreateMatrix4PCF();         
 float CreateMatrix4PCQ();         
 float CreateMatrix4LPCF();         
 float CreateMatrix4LPCQ();         
 float CreateMatrix4PICASSO();         
 float CreateMatrix4PICASSOQ();         
 float CreateMatrix4COMPASS();         
 float CreateMatrix4COMPASSW();         
 float CreateMatrix4PROFSIM(); 
 float CreateMatrix4max(); 
 float      lambda();
 float     *fvector();
 float    **Fmatrix();
 void       free_Fmatrix();
 void       free_fvector();

  

   lena=template[0].main.length;
   lenb=target[0].main.length;
   for (i=0;i<lena;i++){for (j=0;j<lenb;j++){matrix[j][i]=0.0;}}
   for (i=0;i<lena;i++){g[i] =0.0; ge[i]=0.0;}
   for (j=0;j<lenb;j++){p[j] =0.0; pe[j]=0.0;}
   for (k=0;k<MAXP;k++){cte[k]=cte0[k];}
   for (k=0;k<MAXM;k++){weight[k]=weight0[k];}

   maxim=0.0;
          switch (method) {
            case 1:
               maxim=CreateMatrix4raw(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,weight,data,cte,lambda1,tolerance);
               break;
            case 2:
               maxim=CreateMatrix4SF(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,weight,data,cte,lambda1,tolerance);
               break;
            case 3:
               maxim=CreateMatrix4SQ(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
            case 4:
               maxim=CreateMatrix4max(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,weight,data,cte,lambda1,tolerance);
               break;
            case 5:
               maxim=CreateMatrix4PCF(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
            case 6:
               maxim=CreateMatrix4PCQ(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
            case 7:
               maxim=CreateMatrix4LPCF(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
            case 8:
               maxim=CreateMatrix4LPCQ(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
            case 9:
               maxim=CreateMatrix4PICASSO(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
            case 10:
               maxim=CreateMatrix4PICASSOQ(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
            case 11:
               maxim=CreateMatrix4COMPASS(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
            case 12:
               maxim=CreateMatrix4COMPASSW(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
            case 13:
               maxim=CreateMatrix4PROFSIM(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
            default:
               maxim=CreateMatrix4COMPASSW(matrix,g,ge,p,pe,setpro,setmat,checkpro,seq,template,target,PSICtemplate,PSICtarget,qtemplate,qtarget,weight,data,cte,lambda1,tolerance);
               break;
         }

      printf ("Maximum Score in METHOD[%d]: %f\n",method,maxim);
  
      return  maxim;

}
