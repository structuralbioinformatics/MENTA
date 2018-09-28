#include  "inc.h"
/*  PROGRAM  MENTA:  Multiple Entities Alignment */
/*  Baldomero Oliva. Computational Structural Biology Laboratory */
/*  Universitat Pompeu Fabra */
/*  Barcelona. Catalonia, Spain (EU */


void LocalTrace(m,a,b,g,ge,p,pe,table,border,dimension,dimensionj,lms,tolerance,pssm,r)
float      **m;
sequence    a;
sequence    b;
float       *g,*ge,*p,*pe;
tabulation  **table;
boundary    *border;
int         dimension;
int         *dimensionj;
int         lms;
float       *tolerance;
float      **pssm;
alignment   *r;
{

/* work variables */

 int      i,j,k,n,d,loop,dtmp,match,upper,left,h,hh,s;
 int      dj,dd,itmp,jtgt;
 sequence *template,*target;
 reference *Itemplate,*Itarget;
 int      *ii,*jj,*length,*exist;
 float    *score,*penalty;



/* functions */

 sequence    *svector();
 reference   *rvector();
 int         *ivector();
 float       *fvector();
 void         free_svector();
 void         free_rvector();
 void         free_ivector();
 void         free_fvector();


 if (dimension>DIM){printf("Re-dimension Local Subalignments from %d to DIM\n",dimension); dimension=DIM;}

 template=svector(0,dimension+10);
 target=svector(0,dimension+10);
 Itemplate=rvector(0,dimension+10);
 Itarget=rvector(0,dimension+10);
 ii=ivector(0,dimension+10);
 jj=ivector(0,dimension+10);
 length=ivector(0,dimension+10);
 score=fvector(0,dimension+10);
 penalty=fvector(0,dimension+10);
 exist=ivector(0,dimension+10);

 for (i=0;i<dimension;i++){
     exist[i]=0;
 for (j=0;j<MAXS;j++){
     Itemplate[i].image[j]= -1;
     Itarget[i].image[j]= -1;
     Itemplate[i].position[j]= -1;
     Itarget[i].position[j]= -1;
     template[i].element[j]=0;
     target[i].element[j]=0;
 }}
 dtmp=1;
 ii[0]=border[0].imax;
 jj[0]=border[0].jmax;
 length[0]=0;
 for (d=0;d<dimension+10;d++){score[d]=0.0;penalty[d]=0.0;}
 loop=0;
 while (loop<dimension)
 {
   for (d=0;d<dtmp;d++){
    if (dtmp>=dimension+1){loop=dimension;break;}
    i=ii[d];
    j=jj[d];
    k=length[d];
    if ( i<=0 || j<=0) {loop++;break;}
    if ( i>=border[0].imin || j>=border[0].jmin ){
     match=table[j][i].match;
     upper=table[j][i].up;
     left =table[j][i].left;
     if (match == 0 && upper==0 && left==0) {loop++;break;}
     if (exist[d]==0){exist[d]=1;}
     if (match == 1 && upper==0 && left==0) {
       template[d].element[k]=a.element[i];
       target[d].element[k]=b.element[j];
       Itemplate[d].image[i]=j;
       Itemplate[d].position[k]=i;
       Itarget[d].image[j]=i;
       Itarget[d].position[k]=j;
       if (m[j][i]>MIN) score[d]+=m[j][i];
       jj[d]--;
       ii[d]--;
       length[d]++;
     }
     if (match == 0 && upper==1 && left==0) {
       template[d].element[k]=0;
       target[d].element[k]=b.element[j];
       Itarget[d].image[j]= -1;
       Itarget[d].position[k]=j;
       if (penalty[d]==0.0){penalty[d]= -p[j];}else{penalty[d]= -pe[j];}
       if (penalty[d]>MIN) score[d]+= penalty[d];
       jj[d]--;
       length[d]++;
       j=jj[d];
     }
     if (match == 0 && upper==0 && left==1) {
       template[d].element[k]=a.element[i];
       target[d].element[k]=0;
       Itemplate[d].image[i]= -1;
       Itemplate[d].position[k]=i;
       if (penalty[d]==0.0){penalty[d]= -g[i];}else{penalty[d]= -ge[i];}
       if (penalty[d]>MIN) score[d]+= penalty[d];
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
       if (m[j][i]>MIN)score[d]+=m[j][i];
       penalty[d] = 0.0;
       jj[d]--;
       ii[d]--;
       length[d]++;
       template[dtmp].element[k]=0;
       target[dtmp].element[k]=b.element[j];
       Itarget[dtmp].image[j]= -1;
       Itarget[dtmp].position[k]=j;
       if (penalty[dtmp]==0.0){penalty[dtmp]= -p[j];}else{penalty[dtmp]= -pe[j];}
       if (penalty[dtmp]>MIN) score[dtmp]+= penalty[dtmp];

       jj[dtmp]--;
       length[dtmp]++;
/* Total of sequences is increased */
       dtmp++;
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
       if (m[j][i]>MIN)score[d]+=m[j][i];
       penalty[d] = 0.0;
       jj[d]--;
       ii[d]--;
       length[d]++;
       template[dtmp].element[k]=a.element[i];
       target[dtmp].element[k]=0;
       Itemplate[dtmp].image[i]= -1;
       Itemplate[dtmp].position[k]=i;
       if (penalty[dtmp]==0.0){penalty[dtmp]= -g[i];}else{penalty[dtmp]= -ge[i];}
       if (penalty[dtmp]>MIN)score[dtmp]+= penalty[dtmp];
       ii[dtmp]--;
       length[dtmp]++;
/* Total of sequences is increased */
       dtmp++;
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
       ii[dtmp]=ii[d];
       jj[dtmp]=jj[d];
       length[dtmp]=length[d];
/* Add new elements on d and dtmp */
       template[d].element[k]=0;
       target[d].element[k]=b.element[j];
       Itarget[d].image[j]= -1;
       Itarget[d].position[k]=j;
       if (penalty[d]==0.0){penalty[d]= -p[j];}else{penalty[d]= -pe[j];}
       if (penalty[d]>MIN)score[d]+= penalty[d];

       jj[d]--;
       length[d]++;
       template[dtmp].element[k]=a.element[i];
       target[dtmp].element[k]=0;
       Itemplate[dtmp].image[i]= -1;
       Itemplate[dtmp].position[k]=i;
       if (penalty[dtmp]==0.0){penalty[dtmp]= -g[i];}else{penalty[dtmp]= -ge[i];}
       if (penalty[dtmp]>MIN)score[dtmp]+= penalty[dtmp];
       ii[dtmp]--;
       length[dtmp]++;
/* Total of sequences is increased */
       dtmp++;
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
       ii[dtmp]=ii[d];
       jj[dtmp]=jj[d];
       length[dtmp]=length[d];
       score[dtmp+1]=score[d];
       penalty[dtmp+1]=penalty[d];
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
       if (m[j][i]>MIN)score[d]+=m[j][i];
       penalty[d] = 0.0;
       jj[d]--;
       ii[d]--;
       length[d]++;
/* DTMP */
       template[dtmp].element[k]=0;
       target[dtmp].element[k]=b.element[j];
       Itarget[dtmp].image[j]= -1;
       Itarget[dtmp].position[k]=j;
       if (penalty[dtmp]==0.0){penalty[dtmp]= -p[j];}else{penalty[dtmp]= -pe[j];}
       if (penalty[dtmp]>MIN)score[dtmp]+= penalty[dtmp];
       jj[dtmp]--;
       length[dtmp]++;
/* DTMP+1 */
       template[dtmp+1].element[k]=a.element[i];
       target[dtmp+1].element[k]=0;
       Itemplate[dtmp+1].image[i]= -1;
       Itemplate[dtmp+1].position[k]=i;
       if (penalty[dtmp+1]==0.0){penalty[dtmp+1]= -g[i];}else{penalty[dtmp+1]= -ge[i];}
       score[dtmp+1]+= penalty[dtmp+1];
       ii[dtmp+1]--;
       length[dtmp+1]++;
/* Total of sequences is increased twice */
       dtmp=dtmp+2;
      }
    }else{loop++;}
   }
 }
 dd=0;

 *dimensionj=dtmp;


 for (dj=0;dj<dimension;dj++){
   if (dj >= *dimensionj) break; 
   d=dj+dd;
   if (exist[d]==0)break;
   k=length[d];
   h=0;
   hh=0;
   while ( template[d].element[h] == 0
          || target[d].element[h] == 0  ) {h++; if (h>MAXS)break;}
   while ( template[d].element[k-hh] == 0
          || target[d].element[k-hh] == 0  ) {hh++; if (k-hh < 0 ) break;}
   if (  ( (k-h-hh+1) > lms && (score[d] - (hh+h)*(g[0]+p[0])/2)>tolerance[3] && (k-h-hh+1)<MAXS )
       ||( (k-h-hh+1) > lms && tolerance[3] < -1.0e+50 && (k-h-hh+1)<MAXS ) ) {
    r[dj].score=score[d] - (hh+h)*(g[0]+p[0])/2;
    r[dj].evalue=0.0;
    r[dj].pvalue=0.0;
    memset(r[dj].title,'\0',MAXTITLE);
    sprintf(r[dj].title,"ALIGNMENT %5d SCORE %.3e \n",d,r[d].score);
    strcpy(r[dj].template.title,a.title);
    strcpy(r[dj].target.title,b.title);
    i=j=s=0;
    while (j<k-h-hh+1){
     n=k-j-hh;
     if (!(Itemplate[d].position[n] < 0 && Itarget[d].position[n] < 0) ){
      r[dj].template.element[i]   =template[d].element[n];
      r[dj].target.element[i]     =target[d].element[n];
      r[dj].Itemplate.position[i] =Itemplate[d].position[n];
      r[dj].Itarget.position[i]   =Itarget[d].position[n];
      itmp =Itemplate[d].position[n];
      jtgt =Itarget[d].position[n];
      r[dj].Itemplate.image[itmp] =Itemplate[d].image[itmp];
      r[dj].Itarget.image[jtgt]   =Itarget[d].image[jtgt];
      r[dj].profile[i].feature[0] ='\0';
      r[dj].profile[i].feature[1] ='\0';
      r[dj].profile[i].feature[2] ='\0';
      r[dj].profile[i].feature[3] ='\0';
      r[dj].profile[i].feature[4] ='\0';
      if (target[d].element[n]>0 && template[d].element[n]>0 ) pssm[target[d].element[n]][template[d].element[n]] += 1.0;
      i++;
     }else{
      s++;
     }
     j++;
    }
    r[dj].template.length=k-h-hh+1-s;
    r[dj].target.length=k-h-hh+1-s;
    r[dj].length=k-h-hh+1-s;
   }else{
    dd++;
    *dimensionj= dtmp - dd;
    dj--;
   }
 }
 free_ivector(ii,0,dimension+10);
 free_ivector(jj,0,dimension+10);
 free_ivector(exist,0,dimension+10);
 free_ivector(length,0,dimension+10);
 free_svector(template,0,dimension+10);
 free_svector(target,0,dimension+10);
 free_rvector(Itemplate,0,dimension+10);
 free_rvector(Itarget,0,dimension+10);
 free_fvector(score,0,dimension+10);
 free_fvector(penalty,0,dimension+10);


}
