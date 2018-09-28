#include  "inc.h"
/*  PROGRAM  MENTA:  Multiple Entities Alignment */
/*  Baldomero Oliva. Structural Bioinformatics Laboratory */
/*  Universitat Pompeu Fabra */
/*  Barcelona. Catalonia, Spain (EU) */

int   EVAL=0;

main(int argc, char *argv[])
{

/* Variables */
 FILE     *OUT;
 char     out[MAXS],title[MAXP][MAXQ][MAXTITLE],seq[MAXP][MAXQ][MAXS],file[MAXM][500],traduce[MAXE];
 int   i,j,k,n,nn,lms,method,format,psic_flag,lena,lenb,ii,jj,setseq,setmat,setpro,checkpro,dimension,number_of_profiles,old_profiles,skip,global,local,help,*n_set,*id_set,**set;
 float tolerance[MAXTOL],tolerance_pair[MAXTOL],**matrix,**pssm,**pssm_method,*g,*ge,*p,*pe,data[2][MAXM][MAXE][MAXE],maxim,maximum,**qtemplate,**qtarget;
 float weight[MAXM],weight_method[MAXMETH],cte[MAXP],weight_pair[MAXM],cte_pair[MAXP];
 profile  *template,*target,*tgtmp;
 sequence Stemplate,Starget;
 psic   **PSICtemplate,**PSICtarget;
 alignment  *result,dummy_align;
 evd        evd_param;
/* Functions */
 char      *cvector();
 int      LocalDyna(),GlobalDyna(),ReadInput(),OrderAlign(),*ivector(),**imatrix();
 float    CreateMatrix4PAIR(),CreateMatrix4METHOD(),**Fmatrix(),*fvector(),PSSM2Matrix();
 psic   **MEntPSIC();
 alignment *avector();
 sequence  *svector();
 profile   *pvector(),**pmatrix(),**profiles;
 evd     EVDmenta();
 void    Assign(),WriteOutput(),Align4Profile(),SelectProfiles(),WriteProfile(),WriteProfileK();
 void    free_Fmatrix(),free_avector(),free_cvector(),free_svector(),free_fvector(),free_pvector(),free_ivector(),free_imatrix(),free_PSICmatrix(),free_pmatrix();
 void    exit();
 

  setmat=setpro=setseq=checkpro=0;

  id_set=ivector(0,MAXF);
  n_set=ivector(0,MAXF);
  set=imatrix(0,MAXF,0,MAXQ);
  global=local=help=0;
  
  for (i=0;i<argc;i++){
     if (strcmp(argv[i], "-g") == 0) {global=1;}	/*Menta Run in GLOBAL mode  */
  }
  for (i=0;i<argc;i++){
     if (strcmp(argv[i], "-l") == 0) {local=1;}		/*Menta Run in LOCAL mode */
  }
  for (i=0;i<argc;i++){
     if (strcmp(argv[i], "-h") == 0) {help=1;}		/*Menta HELP */
  }
  if (help==1){
    printf("\nPARAMETERS\n");
    printf("\
    \n\t -i \t File with Sequences  \
    \n\t -o \t Output File \
    \n\t -op\t Minimum number of sequences on profiles to write in the output \
    \n\t -a \t Elements (the first two for gap/extension)<default:--ABCDEFGHIKLMNPQRSTVWXYZ>  \
    \n\t -w \t File with Weights and Substitution-Matrix files  \
    \n\t -wj\t File with Weights for Scoring Method (1..13 weights) \
    \n\t -n \t Minimum Alignment <default: 2> \
    \n\t -r \t Tolerance Error <default: 1.0e-6> \
    \n\t -e \t Tolerance Extension  <default: 10 positions> \
    \n\t -ej\t Tolerance Extension for local alignment \
    \n\t    \t only if Scoring-Method \"0\" is applied  <default: NULL => 500 positions > \
    \n\t -ep\t Tolerance Extension in score penalty  <default: 100.0 > \
    \n\t -ejp\t Tolerance Extension penalty for local alignment \
    \n\t    \t only if Scoring-Method \"0\" is applied  <default: NULL => 500*median-score> \
    \n\t -f \t Tolerance Factor <default: 0.0> \
    \n\t -s \t Tolerance Score <default: -100.0> \
    \n\t -id\t Minimum Ratio of sequence identity to create a new profile \
    \n\t -homo\t Minimum Ratio of sequence similarity to create a new profile \
    \n\t -cluster\t Maximum number of profiles with common sequences \
    \n\t -gid\t Gradient factor to apply on restrictions of ID ratio in each iteration\
    \n\t -ghom\t Gradient factor to apply on restrictions of HOMO ratio in each iteration\
    \n\t -psic\t Pseudo Counting flag: ON (simplified approach) / OFF (full pseudo-counting) <default OFF> \
    \n\t -sub\t Level of Global subalignments <default: 3> \
    \n\t    \t   0. Get all subalignments  \
    \n\t    \t   1. Switch off additional pairs of gaps in the same position \
    \n\t    \t   2. Switch off additional gaps in one direction \
    \n\t    \t   3. Switch off all additional gaps (discard subalignments) \
    \n\t -j \t Scoring Method <default: 0> \
    \n\t    \t   0. MEntA \
    \n\t    \t   1. Sum of pairs with Raw Frequencies \
    \n\t    \t   2. Sum of pairs with Efficient Frequencies (F)  \
    \n\t    \t   3. Sum of pairs with Target Frequencies (Q)  \
    \n\t    \t   4. Maximum pair  \
    \n\t    \t   5. Pearson Correlation of Efficient Frequencies (F) \
    \n\t    \t   6. Pearson Correlation of Target Frequencies (Q) \
    \n\t    \t   7. Pearson Correlation of Natural logarithm of F/q \
    \n\t    \t   8. Pearson Correlation of Natural logarithm of Q/q \
    \n\t    \t   9. PICASSO \
    \n\t    \t  10. PICASSO-Q \
    \n\t    \t  11. COMPASS (unweighted) \
    \n\t    \t  12. COMPASS (weighted) \
    \n\t    \t  13. PROFSIM  \
    \n\t -fmt \t Output Format <default: 0> \
    \n\t    \t   0. MEntA format \
    \n\t    \t   1. PIR format   \
    \n\t    \t   2. CLUSTAL format  \
    \n\t    \t  10. MEntA format (includes the alignment of properties) \
    \n\t    \t  11. PIR format (includes the alignment of properties) \
    \n\t    \t  12. CLUSTAL format (includes the alignment of properties) \
    \n\t -evd \t Method of E-value calculation <default: 0> \
    \n\t    \t   0. MEntA  \
    \n\t    \t   1. Island \
    \n\t    \t   2. Murzin P-value \
    \n\t    \t   3. Dummy \
    \n\t -g \t Global Mode  \
    \n\t -l \t Local Mode  \
    \n\t -h \t Print help \n\n");
  }
  if (local==0 && global==0){ printf("Choose running mode (Global= -g | Local= -l | Help= -h)\n\n"); exit(0);}

/* Output file */
  memset(out,'\0',MAXS);
  for (i=0;i<argc;i++){ if (strcmp(argv[i], "-o") == 0) {strcpy(out,argv[i+1]);}}
  if (strcmp(out,"")){OUT=fopen(out,"w");}else{OUT=stdout;}

/* Input Data */
  setseq=ReadInput(argc,argv,&lms,&method,&format,&psic_flag,tolerance,traduce,&setmat,&setpro,title,seq,file,data,weight,weight_method,cte,&checkpro,n_set,set,id_set,&dimension);

/* Allocate Profiles */
  template=pvector(0,setpro);
  target  =pvector(0,setpro);
  tgtmp   =pvector(0,setpro);
  profiles=pmatrix(0,MAXPROF,0,setpro);
  number_of_profiles=0; 

/* Generate profiles */
  for (ii=0;ii<dimension;ii++){
      Assign(seq,title,ii,traduce,setpro,n_set,id_set,set,template,number_of_profiles);
      for (j=0;j<setpro;j++){profiles[ii][j]=template[j];}
      number_of_profiles++;
  }

/*  Print DATA   */
  
  printf("\nSequences are:\n");
  for (i=0;i<setseq;i++){
    printf("\n#%2d > (%d) %s \n",i,strlen(seq[0][i]),title[0][i]);
    printf("     %s\n",seq[0][i]);
  }
  if (setpro>1) printf("\nProperties are:\n");
  for (j=1;j<setpro;j++){
    printf("\n#%2d > Property \n",j);
   for (i=0;i<setseq;i++){
    printf(" $%2d : Score: %e %s \n",i,cte[j-1],title[j][0]);
    printf("     %s\n",seq[j][i]);
   }
  }

  fprintf(OUT,"\n****** START PROFILES  *********\n"); 
  printf("\n****** START PROFILES  *********\n"); 
  for (ii=0;ii<dimension;ii++){
      for (j=0;j<setpro;j++){  WriteProfileK(OUT,traduce,profiles,j,ii); }
  }
  fprintf(OUT,"\n**************************\n"); 
  printf("\n**************************\n"); 
      
/***  ++++++++++++++++++    ***/

  old_profiles=0;

  while (number_of_profiles>old_profiles){

  dimension=number_of_profiles;
    
   for (ii=0;ii<dimension;ii++){
   for (jj=ii+1;jj<dimension;jj++){
    skip=0;
    if (profiles[ii][0].action==0 || profiles[jj][0].action==0) break;
    for (k=0;k<profiles[ii][0].cluster_dimension;k++){ 
    for (j=0;j<profiles[jj][0].cluster_dimension;j++){
      if (profiles[ii][0].cluster[k]==profiles[jj][0].cluster[j]){skip=1;break;}
    }}
    if (jj<old_profiles){skip=1;}
    if (skip==0){

    

     for (j=0;j<setpro;j++){template[j]=profiles[ii][j];}
     for (j=0;j<setpro;j++){target[j]  =profiles[jj][j];}
     lena=template[0].main.length;
     lenb=target[0].main.length;
     Stemplate=template[0].main;
     Starget  =target[0].main;

     matrix=Fmatrix(0,lenb,0,lena);
     pssm  =Fmatrix(0,lenb+1,0,lena+1);
     g     =fvector(0,lena);
     ge    =fvector(0,lena);
     p     =fvector(0,lenb);
     pe    =fvector(0,lenb);

     n=nn=0;

     
/*
     printf("\n> TEMPLATE PROFILE\n");
     fprintf(OUT,"\n> TEMPLATE PROFILE\n");
     for (j=0;j<setpro;j++){  
          WriteProfileK(OUT,traduce,template,j);
       }
     printf("\n**************************\n"); 
     printf("\n> TARGET PROFILE\n");
     fprintf(OUT,"\n**************************\n"); 
     fprintf(OUT,"\n> TARGET PROFILE\n");
     for (j=0;j<setpro;j++){  WriteProfileK(OUT,traduce,target,j); }
     fprintf(OUT,"\n**************************\n"); 
     for (j=0;j<setpro;j++){  WriteProfileK(OUT,traduce,profiles[jj],j); }
     fprintf(OUT,"\n**************************\n"); 

*/


      for(k=0;k<MAXTOL;k++)tolerance_pair[k]=tolerance[k]; /* refresh values */
      for(k=0;k<MAXP;k++)cte_pair[k]=cte[k];               /* refresh values */
      for(k=0;k<MAXM;k++)weight_pair[k]=weight[k];         /* refresh values */

      maxim=CreateMatrix4PAIR(pssm,matrix,g,ge,p,pe,lms,method,psic_flag,tolerance_pair,setpro,setmat,checkpro,seq,template,target,weight_pair,weight_method,data,cte_pair);
      printf("\n**************************\n"); 
      maximum= (lms*maxim)<=MAXV?(lms*maxim):MAXV;
      for(k=0;k<MAXTOL;k++)tolerance_pair[k]=tolerance[k];  /* Use correct tolerances for EVD */
      evd_param=EVDmenta(matrix,Stemplate,Starget,g,ge,p,pe,lms,maximum,tolerance_pair);
      result=avector(0,DIM);
      if (local==1){
        for(k=0;k<MAXTOL; k++)tolerance_pair[k]=tolerance[k]; /* Use correct Tolerance for Local Alignment */
        maximum= (lms*maxim)<=MAXV?(lms*maxim):MAXV;
        printf("\n**************************\n"); 
        printf("Highest Score[Max=%e]: %e\n",maxim,maximum); 
        n=LocalDyna(matrix,Stemplate,Starget,g,ge,p,pe,lms,evd_param,maximum,tolerance_pair,pssm,result);
        printf("Local alignments: %d\n",n);
        if (n>0){nn=OrderAlign(method,traduce,template,target,n,matrix,result,tolerance);}else{nn=0;}
        printf("Local alignments[%d][%d] in profiles: %d\n",ii,jj,nn);
      }
      if (global==1){
        for(k=0;k<MAXTOL; k++)tolerance_pair[k]=tolerance[k]; /* Use correct Tolerance for Global Alignment */
        maximum= (lena*maxim)<=(lenb*maxim)?(lena*maxim):(lenb*maxim);
        printf("\n**************************\n"); 
        printf("Highest Score[Max=%e]: %e\n",maxim,maximum); 
        n=GlobalDyna(matrix,Stemplate,Starget,g,ge,p,pe,evd_param,maximum,tolerance_pair,pssm,result);
        printf("Global alignments: %d\n",n);
        if (n>0){nn=OrderAlign(method,traduce,template,target,n,matrix,result,tolerance);}else{nn=0;}
        printf("Global alignments[%d][%d] in profiles: %d\n",ii,jj,nn);
      }
      WriteOutput(OUT,setpro,method,format,traduce,template,target,n,matrix,result);
      for (k=0;k<nn;k++){
       dummy_align=result[k];
       Align4Profile(dummy_align,template,target,tgtmp,matrix,traduce,setpro,number_of_profiles);
       for (j=0;j<setpro;j++){profiles[number_of_profiles][j]=tgtmp[j];}
       printf("PROFILE[%d] constructed with %s\n",number_of_profiles,dummy_align.title);
       number_of_profiles++;
      }
      free_avector(result,0,DIM);


    free_Fmatrix(matrix,0,lenb,0,lena);
    free_Fmatrix(pssm,0,lenb+1,0,lena+1);
    free_fvector(g,0,lena);
    free_fvector(ge,0,lena);
    free_fvector(p,0,lenb);
    free_fvector(pe,0,lenb);

   }}}
   old_profiles=dimension; 
   SelectProfiles(dimension,number_of_profiles,setpro,tolerance,profiles);
   //printf("Number of Profiles: %d => %d \n",old_profiles,number_of_profiles);
   tolerance[8]=tolerance[12]*tolerance[8];
   if (tolerance[8]>0 && tolerance[8]<MINP)tolerance[8]=0.0;
   tolerance[9]=tolerance[13]*tolerance[9];
   if (tolerance[9]>0 && tolerance[8]<MINP)tolerance[9]=0.0;
  }

  printf("\n****   PROFILES  *****\n");
  fprintf(OUT,"\n****   PROFILES  *****\n");

  for (k=0;k<number_of_profiles;k++){
    if (profiles[k][0].action==1 && profiles[k][0].size>=tolerance[11]){
      for (j=0;j<setpro;j++){  WriteProfileK(OUT,traduce,profiles,j,k); }
    }
  }

    close(OUT);
    free_pmatrix(profiles,0,MAXPROF,0,setpro);
    free_ivector(n_set,0,MAXF);
    free_ivector(id_set,0,MAXF);
    free_imatrix(set,0,MAXF,0,MAXQ);
    free_pvector(template,0,setpro);
    free_pvector(target,0,setpro);
    free_pvector(tgtmp,0,setpro);
}
