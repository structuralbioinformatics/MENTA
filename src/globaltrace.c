#include "inc.h"
/* Subroutine developed by B. Oliva  for MEnTA. */

int GlobalTrace(m,a,b,g,ge,p,pe,tolerance,table,dimension,pssm,r)
float      **m;
sequence    a;
sequence    b;
float       *g,*ge,*p,*pe,*tolerance;
tabulation  **table;
int         dimension;
float       **pssm;
alignment   *r;
{


/* work variables */

 int      i,j,k,n,d,loop,dtmp,match,upper,left,h,hh;
 int      itmp,jtgt; 
 int      dd,dimension_2;
 sequence *template,*target;
 reference *Itemplate,*Itarget;
 int      *ii,*jj,*length;
 float    *score,*penalty,*extension;

/* functions */

 sequence    *svector();
 int         *ivector();
 reference   *rvector();
 float       *fvector();
 void         free_svector();
 void         free_rvector();
 void         free_ivector();
 void         free_fvector();

 if (dimension>DIM || dimension<0) return 0;
 dimension_2=dimension+2;

 printf ("Global-Trace Dimension: %d \n",dimension); 
 template=svector(0,dimension_2);
 target=svector(0,dimension_2);
 Itemplate=rvector(0,dimension_2);
 Itarget=rvector(0,dimension_2);
 ii=ivector(0,dimension_2);
 jj=ivector(0,dimension_2);
 length=ivector(0,dimension_2);
 score=fvector(0,dimension_2);
 penalty=fvector(0,dimension_2);
 extension=fvector(0,dimension_2);
 for (i=0;i<dimension_2;i++){
 for (j=0;j<MAXS;j++){
     Itemplate[i].image[j]= -1;
     Itarget[i].image[j]= -1;
     Itemplate[i].position[j]= -1;
     Itarget[i].position[j]= -1;
 }}
 dtmp=1;
 ii[0]=a.length-1;
 jj[0]=b.length-1;
 length[0]=0;
 for (d=0;d<dimension_2;d++){score[d]=0.0;penalty[d]=0.0;extension[d]=0.0;}
 loop=0;
 while (loop<dimension)
 {
   for (d=0;d<dtmp;d++){

    i=ii[d];
    j=jj[d];
    k=length[d];


    if (i>=0 || j>=0){
     match=table[j][i].match;
     upper=table[j][i].up;
     left =table[j][i].left;
     if (i<0){ template[d].element[k]=0; }
     if (j<0){ target[d].element[k]=0; }
     if (match == 1 && upper==0 && left==0) {
       template[d].element[k]=a.element[i];
       target[d].element[k]=b.element[j];
       Itemplate[d].image[i]=j;
       Itemplate[d].position[k]=i;
       Itarget[d].image[j]=i;
       Itarget[d].position[k]=j;
       score[d]+=m[j][i];
       if (extension[d]>=tolerance[2]){score[d]-=extension[d];}
       penalty[d] = 0.0;
       extension[d]=0.0;
       jj[d]--;
       ii[d]--;
       length[d]++;
     }
     if (match == 0 && upper==1 && left==0) {
       template[d].element[k]=0;
       target[d].element[k]=b.element[j];
       Itarget[d].image[j]= -1;
       Itarget[d].position[k]=j;
       if (penalty[d]==0.0 ){ 
	       if (extension[d]<tolerance[2]) {penalty[d]= -p[j];}
       }else{
	       if (extension[d]<tolerance[2]) { penalty[d] = -pe[j];
		                            extension[d] += pe[j];
	       }else{ penalty[d] = 0.0; }
       }
       jj[d]--;
       length[d]++;
       j=jj[d];
     }
     if (match == 0 && upper==0 && left==1) {
       template[d].element[k]=a.element[i];
       target[d].element[k]=0;
       Itemplate[d].image[i]= -1;
       Itemplate[d].position[k]=i;
       if (penalty[d]==0.0){ 
	       if (extension[d]<tolerance[2]) {penalty[d]= -g[i];}
       }else{
	       if (extension[d]<tolerance[2]) { penalty[d] = -ge[i];
		                          extension[d] += ge[i];
	       }else{ penalty[d] = 0.0; }
       }
       ii[d]--;
       length[d]++;
       i=ii[d];
     }
     if (match == 1 && upper==1 && left==0) {
/* Copy sequence of d to dtmp */
       for (n=0;n<length[d];n++){
        template[dtmp].element[n]=template[d].element[n];
        target[dtmp].element[n]=target[d].element[n];
       }
       for (n=0;n<length[d];n++){
        itmp=Itemplate[d].position[n];
        Itemplate[dtmp].position[n]=Itemplate[d].position[n];
        Itemplate[dtmp].image[itmp]=Itemplate[d].image[itmp];
       }
       for (n=0;n<length[d];n++){
        jtgt=Itarget[d].position[n];
        Itarget[dtmp].position[n]=Itarget[d].position[n];
        Itarget[dtmp].image[jtgt]=Itarget[d].image[jtgt];
       }
       score[dtmp]=score[d];
       penalty[dtmp]=penalty[d];
       extension[dtmp]=extension[d];
       ii[dtmp]=ii[d];
       jj[dtmp]=jj[d];
       length[dtmp]=length[d];
/* Add new elements on d and dtmp */
       template[d].element[k]=a.element[i];
       target[d].element[k]=b.element[j];
       Itemplate[d].image[i]=j;
       Itemplate[d].position[k]=i;
       Itarget[d].image[j]=i;
       Itarget[d].position[k]=j;
       score[d]+=m[j][i];
       if (extension[d]>=tolerance[2]){score[d]-=extension[d];}
       penalty[d] = 0.0;
       extension[d]=0.0;
       jj[d]--;
       ii[d]--;
       length[d]++;
       template[dtmp].element[k]=0;
       target[dtmp].element[k]=b.element[j];
       Itarget[dtmp].image[j]= -1;
       Itarget[dtmp].position[k]=j;
       if (penalty[dtmp]==0.0){
	       if (extension[d]<tolerance[2]) {penalty[dtmp]= -p[j];}
       }else{
	       if (extension[dtmp]<tolerance[2]) { penalty[dtmp] = -pe[j];
		                             extension[dtmp] += pe[j];
	       }else{ penalty[dtmp] = 0.0; }
       }
       score[dtmp]+= penalty[dtmp];
       jj[dtmp]--;
       length[dtmp]++;
/* Total of sequences is increased */
       dtmp++;
       if (dtmp>dimension)dtmp=dimension;
     }
     if (match == 1 && upper==0 && left==1) {
/* Copy sequence of d to dtmp */
       for (n=0;n<length[d];n++){
        template[dtmp].element[n]=template[d].element[n];
        target[dtmp].element[n]=target[d].element[n];
       }
       for (n=0;n<length[d];n++){
        itmp=Itemplate[d].position[n];
        Itemplate[dtmp].position[n]=Itemplate[d].position[n];
        Itemplate[dtmp].image[itmp]=Itemplate[d].image[itmp];
       }
       for (n=0;n<length[d];n++){
        jtgt=Itarget[d].position[n];
        Itarget[dtmp].position[n]=Itarget[d].position[n];
        Itarget[dtmp].image[jtgt]=Itarget[d].image[jtgt];
       }
       score[dtmp]=score[d];
       penalty[dtmp]=penalty[d];
       extension[dtmp]=extension[d];
       ii[dtmp]=ii[d];
       jj[dtmp]=jj[d];
       length[dtmp]=length[d];
/* Add new elements on d and dtmp */
       template[d].element[k]=a.element[i];
       target[d].element[k]=b.element[j];
       Itemplate[d].image[i]=j;
       Itemplate[d].position[k]=i;
       Itarget[d].image[j]=i;
       Itarget[d].position[k]=j;
       score[d]+=m[j][i];
       if (extension[d]>=tolerance[2]){score[d]-=extension[d];}
       penalty[d] = 0.0;
       extension[d] = 0.0;
       jj[d]--;
       ii[d]--;
       length[d]++;
       template[dtmp].element[k]=a.element[i];
       target[dtmp].element[k]=0;
       Itemplate[dtmp].image[i]= -1;
       Itemplate[dtmp].position[k]=i;
       if (penalty[dtmp]==0.0){
	       if (extension[d]<tolerance[2]) {penalty[dtmp]= -g[i];}
       }else{
	       if (extension[dtmp]<tolerance[2]) { penalty[dtmp] = -ge[i];
		                             extension[dtmp] += ge[i];
	       }else{penalty[dtmp] = 0.0; }
       }
       score[dtmp]+= penalty[dtmp];
       ii[dtmp]--;
       length[dtmp]++;
/* Total of sequences is increased */
       dtmp++;
       if (dtmp>dimension)dtmp=dimension;
     }
     if (match == 0 && upper==1 && left==1) {
/* Copy sequence of d to dtmp */
       for (n=0;n<length[d];n++){
        template[dtmp].element[n]=template[d].element[n];
        target[dtmp].element[n]=target[d].element[n];
       }
       for (n=0;n<length[d];n++){
        itmp=Itemplate[d].position[n];
        Itemplate[dtmp].position[n]=Itemplate[d].position[n];
        Itemplate[dtmp].image[itmp]=Itemplate[d].image[itmp];
       }
       for (n=0;n<length[d];n++){
        jtgt=Itarget[d].position[n];
        Itarget[dtmp].position[n]=Itarget[d].position[n];
        Itarget[dtmp].image[jtgt]=Itarget[d].image[jtgt];
       }
       score[dtmp]=score[d];
       penalty[dtmp]=penalty[d];
       extension[dtmp]=extension[d];
       ii[dtmp]=ii[d];
       jj[dtmp]=jj[d];
       length[dtmp]=length[d];
/* Add new elements on d and dtmp */
       template[d].element[k]=0;
       target[d].element[k]=b.element[j];
       Itarget[d].image[j]= -1;
       Itarget[d].position[k]=j;
       if (penalty[d]==0.0){
	       if (extension[d]<tolerance[2]) {penalty[d]= -p[j];}
       }else{
	       if (extension[d]<tolerance[2]) {penalty[d] = -pe[j];
		                         extension[d]+=pe[j];
	       }else{penalty[d] = 0; }
       }
       jj[d]--;
       length[d]++;
       template[dtmp].element[k]=a.element[i];
       target[dtmp].element[k]=0;
       Itemplate[dtmp].image[i]= -1;
       Itemplate[dtmp].position[k]=i;
       if (penalty[dtmp]==0.0){
	       if (extension[d]<tolerance[2]) {penalty[dtmp]= -g[i];}
       }else{
	       if (extension[dtmp]<tolerance[2]) {penalty[dtmp] = -ge[i];
		                              extension[dtmp] += ge[i];
	       }else{penalty[dtmp] = 0; }
       }
       ii[dtmp]--;
       length[dtmp]++;
/* Total of sequences is increased */
       dtmp++;
       if (dtmp>dimension)dtmp=dimension;
     }
     if (match == 1 && upper==1 && left==1) {
/* Copy twice sequence of d to dtmp and dtmp+1*/
       for (n=0;n<length[d];n++){
        template[dtmp].element[n]=template[d].element[n];
        target[dtmp].element[n]=target[d].element[n];
        template[dtmp+1].element[n]=template[d].element[n];
        target[dtmp+1].element[n]=target[d].element[n];
       }
       for (n=0;n<length[d];n++){
        itmp=Itemplate[d].position[n];
        Itemplate[dtmp].position[n]=Itemplate[d].position[n];
        Itemplate[dtmp].image[itmp]=Itemplate[d].image[itmp];
        Itemplate[dtmp+1].position[n]=Itemplate[d].position[n];
        Itemplate[dtmp+1].image[itmp]=Itemplate[d].image[itmp];
       }
       for (n=0;n<length[d];n++){
        jtgt=Itarget[d].position[n];
        Itarget[dtmp].position[n]=Itarget[d].position[n];
        Itarget[dtmp].image[jtgt]=Itarget[d].image[jtgt];
        Itarget[dtmp+1].position[n]=Itarget[d].position[n];
        Itarget[dtmp+1].image[jtgt]=Itarget[d].image[jtgt];
       }
       score[dtmp]=score[d];
       penalty[dtmp]=penalty[d];
       extension[dtmp]=extension[d];
       ii[dtmp]=ii[d];
       jj[dtmp]=jj[d];
       length[dtmp]=length[d];
       score[dtmp+1]=score[d];
       penalty[dtmp+1]=penalty[d];
       extension[dtmp+1]=extension[d];
       ii[dtmp+1]=ii[d];
       jj[dtmp+1]=jj[d];
       length[dtmp+1]=length[d];
/* Add new elements on d , dtmp and dtmp+1 */
/* D */
       template[d].element[k]=a.element[i];
       target[d].element[k]=b.element[j];
       Itemplate[d].image[i]=j;
       Itemplate[d].position[k]=i;
       Itarget[d].image[j]=i;
       Itarget[d].position[k]=j;
       score[d]+=m[j][i];
       if (extension[d]>=tolerance[2]){score[d]-=extension[d];}
       penalty[d] = 0.0;
       extension[d] = 0.0;
       jj[d]--;
       ii[d]--;
       length[d]++;
/* DTMP */
       template[dtmp].element[k]=0;
       target[dtmp].element[k]=b.element[j];
       Itarget[dtmp].image[j]= -1;
       Itarget[dtmp].position[k]=j;
       if (penalty[dtmp]==0.0){
	       if (extension[d]<tolerance[2]) {penalty[dtmp]= -p[j];}
       }else{
	       if (extension[dtmp]<tolerance[2]) {penalty[dtmp] = -pe[j];
		                            extension[dtmp] += pe[j];
	       }else{penalty[dtmp] = 0; }
       }
       score[dtmp]+= penalty[dtmp];
       jj[dtmp]--;
       length[dtmp]++;
/* DTMP+1 */
       template[dtmp+1].element[k]=a.element[i];
       target[dtmp+1].element[k]=0;
       Itemplate[dtmp+1].image[i]= -1;
       Itemplate[dtmp+1].position[k]=i;
       if (penalty[dtmp+1]==0.0){
	       if (extension[d]<tolerance[2]) {penalty[dtmp+1]= -g[i];}
       }else{
	       if (extension[dtmp+1]<tolerance[2]) {penalty[dtmp+1] = -ge[i];
		                              extension[dtmp+1]+= ge[i];
	       }else{penalty[dtmp+1] = 0.0;  }
       }
       score[dtmp+1]+= penalty[dtmp+1];
       ii[dtmp+1]--;
       length[dtmp+1]++;
/* Total of sequences is increased twice */
       dtmp=dtmp+2;
       if (dtmp>dimension)dtmp=dimension;
      }
    }else{loop++;}
   }

   for (d=0;d<dtmp;d++){
    i=ii[d];
    j=jj[d];
    k=length[d];
    if (i<0){ template[d].element[k]=0; }
    if (j<0){ target[d].element[k]=0; }
    if (i>0 || j>0){loop--;}
   }

 }

 dd=0; 
 for (d=0;d<dimension;d++){
   h=0;
   hh=0;
   k=length[d];
   while ( template[d].element[h]  <= 0
          && target[d].element[h]  <= 0  ) {h++;}
   while ( template[d].element[k-hh]  <= 0
          && target[d].element[k-hh]  <= 0  ) {hh++;}
   if ((k-h-hh+1)>0 && (k-h-hh+1)<MAXS){
    r[dd].template.length=k-h-hh+1;
    r[dd].target.length=k-h-hh+1;
    r[dd].length=k-h-hh+1;
    r[dd].score= score[d] - (hh+h)*(g[0]+p[0])/2;
    r[dd].evalue=0.0; 
    r[dd].pvalue=0.0; 
    memset(r[dd].title,'\0',MAXTITLE);
    sprintf(r[dd].title,"ALIGNMENT %5d SCORE %.3e \n",d,r[dd].score);
    strcpy(r[dd].template.title,a.title);
    strcpy(r[dd].target.title,b.title);
    r[dd].homol=0.0;
    r[dd].ident=0.0;
    for (i=0;i<k-h-hh+1;i++){
     n=k-i-hh;
     r[dd].template.element[i]   =template[d].element[n];
     r[dd].target.element[i]     =target[d].element[n];
     r[dd].Itemplate.position[i] =Itemplate[d].position[n]-1;
     r[dd].Itarget.position[i]   =Itarget[d].position[n]-1;
     itmp =Itemplate[d].position[n];
     jtgt =Itarget[d].position[n];
     r[dd].Itemplate.image[itmp-1] =Itemplate[d].image[itmp]-1;
     r[dd].Itarget.image[jtgt-1]   =Itarget[d].image[jtgt]-1;
     r[dd].profile[i].feature[0] ='\0';
     r[dd].profile[i].feature[1] ='\0';
     r[dd].profile[i].feature[2] ='\0';
     r[dd].profile[i].feature[3] ='\0';
     r[dd].profile[i].feature[4] ='\0';
     pssm[target[d].element[n]][template[d].element[n]] += 1.0;
    }
    dd++;
   }else{
    printf("Check MAXS sequence length %d > MAXS\n",k-h-hh+1);
   }
 }


 free_ivector(ii,0,dimension_2);
 free_ivector(jj,0,dimension_2);
 free_ivector(length,0,dimension_2);
 free_svector(template,0,dimension_2);
 free_svector(target,0,dimension_2);
 free_rvector(Itemplate,0,dimension_2);
 free_rvector(Itarget,0,dimension_2);
 free_fvector(score,0,dimension_2);
 free_fvector(penalty,0,dimension_2);
 free_fvector(extension,0,dimension_2);
 return dd;
}
