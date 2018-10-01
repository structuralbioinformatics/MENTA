#include  "inc.h"
int ReadInput(argc,argv,lms,method,format,psic_flag,tolerance,traduce,setmat,setpro,title,seq,file,data,weight,weight_method,cte,checkpro,n_set,set,id_set,dimension)
int     argc;
char   *argv[];
int    *lms,*method,*format,*psic_flag;
float  *tolerance;
char   traduce[MAXE];
int    *setmat;
int    *setpro;
char   title[MAXP][MAXQ][MAXTITLE];
char   seq[MAXP][MAXQ][MAXS];
char   file[MAXM][500];
float  data[2][MAXM][MAXE][MAXE];
float  weight[MAXM];
float  weight_method[MAXMETH];
float  cte[MAXP];
int    *checkpro,*dimension;
int    *n_set,**set,*id_set;
{
 
 FILE     *INP;
 char     input[MAXS],inputw[MAXS],inputwj[MAXS],buffer[MAXS],data_traduce[MAXE],check[1];
 char     title_mat[MAXP][MAXQ][MAXTITLE],error_size[MAXS];
 int      np,ns,ii,i,j,k,length[MAXQ],setseq;
 int      id_profile,dim;
 int      tmp_profile[MAXQ],tmp_ii_profile[MAXQ],tmp_check,ii_tmp,ii_profile,got;
 void     nrerror(); 

  memset(input,'\0',MAXS);
  memset(inputw,'\0',MAXS);
  memset(data_traduce,'\0',MAXE);
  memset(check,'\0',1);
  strcpy(data_traduce,"--ABCDEFGHIKLMNPQRSTVWXYZ");
  *lms=2;
  *method=0;
  *format=0;
  *psic_flag=0;

  tolerance[0]=  0.0;
  tolerance[1]=  1.0e-3;
  tolerance[2]=  1.0e+1;
  tolerance[3]= -1.0e+2;
  tolerance[4]=  3.0;
  tolerance[5]=  0.0;
  tolerance[6]=  100.0;
  for (i=7;i<MAXTOL;i++){ tolerance[i]=  0.0;}
  tolerance[10]= 3.0;
  tolerance[11]= 2.0;
  tolerance[12]= 1.0;
  tolerance[13]= 1.0;

  for (i=0;i<MAXP;i++){cte[i]=0.0;}
/* Reading input line */
 
  for (i=0;i<argc;i++){
      if (strcmp(argv[i], "-i") == 0) {	/*File with Sequences  */
         strcpy(input, argv[i+1]);
      }
      if (strcmp(argv[i], "-a") == 0) {	/*Elements (the first two for gap/extension)  */
         strcpy(data_traduce, argv[i+1]);
      }
      if (strcmp(argv[i], "-w") == 0) {	/*File with Weights and Matrix files */
         strcpy(inputw, argv[i+1]);
      }
      if (strcmp(argv[i], "-wj") == 0) {/*File with Weights for Scoring Method */
         strcpy(inputwj, argv[i+1]);
      }
      if (strcmp(argv[i], "-n") == 0) {	/*Minimum Alignment Length  */
         sscanf(argv[i+1], "%d", lms);
      }
      if (strcmp(argv[i], "-f") == 0) {	/*Tolerance Factor  */
         sscanf(argv[i+1], "%f", &tolerance[0]);
      }
      if (strcmp(argv[i], "-r") == 0) {	/*Tolerance error */
         sscanf(argv[i+1], "%f", &tolerance[1]);
      }
      if (strcmp(argv[i], "-e") == 0) {	/*Tolerance Extension */
         sscanf(argv[i+1], "%f", &tolerance[2]);
      }
      if (strcmp(argv[i], "-s") == 0) {	/*Tolerance Score */
         sscanf(argv[i+1], "%f", &tolerance[3]);
      }
      if (strcmp(argv[i], "-sub") == 0) {	/*Tolerance Level of number of subalignments */
         sscanf(argv[i+1], "%f", &tolerance[4]);
      }
      if (strcmp(argv[i], "-ej") == 0) {	/*Tolerance Extension on Local alignment of combined methods*/
         sscanf(argv[i+1], "%f", &tolerance[5]);
      }
      if (strcmp(argv[i], "-ep") == 0) {	/*Tolerance Extension in score penalty*/
         sscanf(argv[i+1], "%f", &tolerance[6]);
      }
      if (strcmp(argv[i], "-ejp") == 0) {	/*Tolerance Extension in score penalty  on Local alignment of combined methods*/
         sscanf(argv[i+1], "%f", &tolerance[7]);
      }
      if (strcmp(argv[i], "-id") == 0) {	/*Tolerance  in percentage of identity  on results of alignment*/
         sscanf(argv[i+1], "%f", &tolerance[8]);
      }
      if (strcmp(argv[i], "-homo") == 0) {	/*Tolerance  in percentage of homology  on results of alignment*/
         sscanf(argv[i+1], "%f", &tolerance[9]);
      }
      if (strcmp(argv[i], "-cluster") == 0) {	/*Tolerance  Maximum number of Profiles with the same sequences*/
         sscanf(argv[i+1], "%f", &tolerance[10]);
      }
      if (strcmp(argv[i], "-op") == 0) {	/*Tolerance  Minimum number of sequences on profiles to write in the output*/
         sscanf(argv[i+1], "%f", &tolerance[11]);
      }
      if (strcmp(argv[i], "-gid") == 0) {	/*Tolerance   gradient factor for restrictions of ID  percentages*/
         sscanf(argv[i+1], "%f", &tolerance[12]);
      }
      if (strcmp(argv[i], "-ghom") == 0) {	/*Tolerance   gradient factor for restrictions of HOMO percentages*/
         sscanf(argv[i+1], "%f", &tolerance[13]);
      }
      if (strcmp(argv[i], "-j") == 0) {	/*Method of alignment */
         sscanf(argv[i+1], "%d", method);
      }
      if (strcmp(argv[i], "-fmt") == 0) {	/*Format of Output-alignment */
         sscanf(argv[i+1], "%d", format);
      }
      if (strcmp(argv[i], "-evd") == 0) {	/*Method of Calculation of e-values */
         sscanf(argv[i+1], "%d", &EVAL);
      }
      if (strcmp(argv[i], "-psic") == 0) {    /*Pseudo counting flag */
          *psic_flag=1;
      }

  } 

  if (*lms==0){*lms=2;}
  if (*method>MAXMETH-1){*method=0;}

  printf("\nRunning Data:\n");
  printf("File with sequences:%s\n",input);
  printf("Minimum alignment: %d \n",*lms);
  printf("Tolerance Factor: %e \n",tolerance[0]);
  printf("Tolerance Error:%e \n",tolerance[1]);
  printf("Limited extension length: %e \n",tolerance[2]);
  printf("Limited extension penalty: %e \n",tolerance[6]);
  printf("Limited local-method extension length: %e \n",tolerance[5]);
  printf("Limited local-method extension penalty: %e \n",tolerance[7]);
  printf("Limited score: %e \n",tolerance[3]);
  printf("Level of Subalignments: %e \n",tolerance[4]);
  printf("Elements (first two must be gap/extension):%s\n",data_traduce);
  printf("File with Weights and Matrix files:%s\n",inputw);
  printf("File with Weights for Scoring Method:%s\n",inputwj);
  printf("Scoring Method: %d \n",*method);
  printf("E-value Method: %d \n",EVAL);
  printf("Minimum percentage of identities to accept alignements: %f \n",100*tolerance[8]);
  printf("Minimum percentage of homology to accept alignements: %f \n",100*tolerance[9]);
  printf("Maximum number of profiles with the same sequences: %f \n",tolerance[10]);
  printf("Minimum number of sequences on profiles to write in the output: %f \n",tolerance[11]);
  printf("Gradient factor to apply on restrictions of ID percentages: %f \n",tolerance[12]);
  printf("Gradient factor to apply on restrictions of HOMO percentages: %f \n",tolerance[13]);

  INP=fopen(input,"r");
  ii=0;
  np=0;
  dim=0;

  if (!INP){nrerror("No data\n");}
  while(!feof(INP)) {
    memset(buffer,'\0',MAXS);
    fgets(buffer,MAXS,INP);
    sscanf(buffer,"%c",check);
    if (check[0]=='>'){ sscanf(buffer+1,"%c",check);
	                if (!strncmp(check,"M",1)) {
                         sscanf(buffer+2,"%d",&id_profile);
	                 memset(title[0][ii],'\0',MAXTITLE);
                         strcpy(title[0][ii],buffer+1);
                         memset(seq[0][ii],'\0',MAXS);
                         length[ii]=0;
                         tmp_profile[ii]=id_profile;
                         tmp_check=0;
                         ii_profile=0;
                         for (ii_tmp=0;ii_tmp<ii;ii_tmp++){
                           if (tmp_profile[ii]==tmp_profile[ii_tmp]){
                               tmp_check=1;
                               ii_profile=tmp_ii_profile[ii_tmp];
                               break;
                           }
                         }
                         if (tmp_check==0){ii_profile=tmp_ii_profile[ii]=dim;dim=dim+1;id_set[ii_profile]=id_profile;}
			 if (dim>=MAXF) nrerror("Too many profiles, change MAXF");
			 set[ii_profile][n_set[ii_profile]]=ii;
			 n_set[ii_profile]= n_set[ii_profile] + 1;
                         if (n_set[ii_profile] >= MAXFS) {
                            sprintf(error_size,"Too Many Sequences on a Profile, change MAXFS: %d",ii_profile);
                            nrerror(error_size);
                         } 
                         ii++;
                         if (ii>=MAXQ) nrerror("Too many sequences,change MAXQ");
                         np=0;
			}else{
	                 memset(title[0][ii],'\0',MAXTITLE);
                         strcpy(title[0][ii],buffer+1);
                         memset(seq[0][ii],'\0',MAXS);
                         length[ii]=0;
                         ii++;
                         if (ii>=MAXQ) nrerror("Too many sequences, change MAXQ");
                         np=0;
			}
    }else if (check[0]=='#') {
                         np++;
                         memset(title_mat[np][ii],'\0',MAXTITLE);
                         if (ii==1 && np>0 ) {
                          sscanf(buffer+1,"%f",&cte[np-1]);
                          sscanf(buffer+10,"%100s",title_mat[np][ii-1]);
                          if (cte[np-1]!=0 && *checkpro==0) *checkpro=1;
                         }
                         memset(seq[np][ii],'\0',MAXS);
                         if (np>=MAXP)  nrerror("Too many properties, change MAXP");}
    else if (check[0]!='\n'){ 
                        ns=strlen(buffer);
                        if (np==0) {length[ii-1] += ns-1;}
                        if (length[ii-1]>=MAXS) {sprintf(error_size,"Too large sequence, change MAXS > %d",length[ii-1]);nrerror(error_size);}
                        strncat(seq[np][ii-1],buffer,ns-1);}
  }
  close(INP);
  setseq=ii;
  *setpro=np+1;
  for (i=0;i<setseq;i++){
    got=0;
    for (j=0;j<dim;j++){
      if (got==0){ for (k=0;k<n_set[j];k++){ if (i==set[j][k]){got=1;break;} } }else{break;} 
    }
    if (got==0){
      set[dim][0]=i;
      n_set[dim]=1;
      dim=dim+1; 
      if (dim>=MAXF) nrerror("Too many profiles, change MAXF");
    }
  }
 
  *dimension=dim;
 
  if (strncmp(data_traduce,"-",1) || strncmp(data_traduce+1,"-",1) ) nrerror("Wrong GAP definition");

  for (i=0;i<strlen(data_traduce);i++){traduce[i]=data_traduce[i];}

  for (i=strlen(data_traduce);i<MAXE;i++){traduce[i]='*';}

  INP=fopen(inputw,"r");
  ii=0;
  if (!INP){nrerror("No W&M Files data\n");}
  while(!feof(INP)) {
    memset(buffer,'\0',MAXS);
    memset(file[ii],'\0',100);
    fgets(buffer,MAXS,INP);
    sscanf(buffer,"%f",&weight[ii]);
    sscanf(buffer+10,"%100s",file[ii]);
    ii++;
    if (ii>=MAXM) nrerror("Too many Matrices");
  }
  close(INP);
  *setmat=ii-1;

/** Read the Matrices of log-odds **/


  for (k=0;k<*setmat;k++){
      INP=fopen(file[k],"r");
       printf ("open %s \n",file[k]);
      if (!INP){
       nrerror ("No MATRIX data\n");
      }else{
       i=0;
       while(!feof(INP)) {
        memset(buffer,'\0',MAXS);
        fgets(buffer,MAXS,INP);
        if (strncmp(buffer,"#",1)) {
         for (j=0;j<MAXE;j++){sscanf(buffer+10*j,"%f",&data[0][k][i][j]);}
         i++;
         if (i>MAXE) { printf ("Row-Residues: %d \n",i); nrerror("Too many elements for a Matrix");}
        }else{
         printf ("\tREMARK %s \n",buffer);
        }
       }
       close(INP);
      }
  }

  for (k=1;k<*setpro;k++){
      INP=fopen(title_mat[k][0],"r");
      printf ("open %s \n",title_mat[k][0]);
      if (!INP){
       nrerror ("No MATRIX data\n");
      }else{
       i=0;
       while(!feof(INP)) {
        memset(buffer,'\0',MAXS);
        fgets(buffer,MAXS,INP);
        if (strncmp(buffer,"#",1)) {
         for (j=0;j<MAXE;j++){sscanf(buffer+10*j,"%f",&data[1][k-1][i][j]);}
         i++;
         if (i>MAXE) { printf ("Row-Properties: %d \n",i); nrerror("Too many elements for a Matrix");}
        }else{
         printf ("\tREMARK %s \n",buffer);
        }
       }
       close(INP);
      }
  }

/** Read the Weights for Scoring **/


  INP=fopen(inputwj,"r");
  ii=0;
  if (!INP){
   printf("No Weights for Method (all=1)\n");
   for (ii=0;ii<MAXMETH;ii++){weight_method[ii]=1.0;}
  }else{
   while(!feof(INP)) {
    memset(buffer,'\0',MAXS);
    fgets(buffer,MAXS,INP);
    sscanf(buffer,"%f",&weight_method[ii]);
    printf ("Weight-Method[%2d]= %f \n",ii+1,weight_method[ii]);
    ii++;
    if (ii>MAXMETH) nrerror("Too many Methods");
   }
   close(INP);
  }

  return setseq;

}
