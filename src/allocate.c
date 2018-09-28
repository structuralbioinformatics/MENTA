#include  "inc.h"


 void nrerror(error_text)
 char error_text[];
 {
 void exit();
 fprintf(stderr,"Numerical Recipes run-time error...\n");
 fprintf(stderr,"%s\n",error_text);
 fprintf(stderr,"...now exiting to system...\n");
 exit(1);
 }


tabulation  **Tmatrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 tabulation  **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(tabulation **)  malloc((unsigned) (nrh-nrl+1)*sizeof(tabulation *));
 if (!m) nrerror("allocation failure 1 in Tmatrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(tabulation *) malloc((unsigned) (nch-ncl+1)*sizeof(tabulation));
  if (!m[i]) nrerror("allocation failure 2 in Tmatrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

void free_Tmatrix(m,nrl,nrh,ncl,nch)
tabulation **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((tabulation *) (m[i]+ncl));
   free((tabulation *) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((int *) (m[i]+ncl));
   free((int *) (m+nrl));
}


int  **imatrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 int  **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(int **)  malloc((unsigned) (nrh-nrl+1)*sizeof(int *));
 if (!m) nrerror("allocation failure 1 in imatrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
  if (!m[i]) nrerror("allocation failure 2 in imatrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

float  **Fmatrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 float  **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(float **)  malloc((unsigned) (nrh-nrl+1)*sizeof(float *));
 if (!m) nrerror("allocation failure 1 in Fmatrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
  if (!m[i]) nrerror("allocation failure 2 in Fmatrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

void free_Fmatrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((float *) (m[i]+ncl));
   free((float *) (m+nrl));
}

psic  **PSICmatrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=PSIC  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 psic  **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(psic **)  malloc((unsigned) (nrh-nrl+1)*sizeof(psic *));
 if (!m) nrerror("allocation failure 1 in PSICmatrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(psic *) malloc((unsigned) (nch-ncl+1)*sizeof(psic));
  if (!m[i]) nrerror("allocation failure 2 in PSICmatrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

void free_PSICmatrix(m,nrl,nrh,ncl,nch)
psic **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((psic *) (m[i]+ncl));
   free((psic *) (m+nrl));
}



char   **cmatrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 char   **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(char  **)  malloc((unsigned) (nrh-nrl+1)*sizeof(char  *));
 if (!m) nrerror("allocation failure 1 in cmatrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(char  *) malloc((unsigned) (nch-ncl+1)*sizeof(char ));
  if (!m[i]) nrerror("allocation failure 2 in cmatrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

void free_cmatrix(m,nrl,nrh,ncl,nch)
char  **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((char  *) (m[i]+ncl));
   free((char  *) (m+nrl));
}

profile   **pmatrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 profile   **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(profile  **)  malloc((unsigned) (nrh-nrl+1)*sizeof(profile  *));
 if (!m) nrerror("allocation failure 1 in pmatrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(profile  *) malloc((unsigned) (nch-ncl+1)*sizeof(profile ));
  if (!m[i]) nrerror("allocation failure 2 in pmatrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

void free_pmatrix(m,nrl,nrh,ncl,nch)
profile  **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((profile  *) (m[i]+ncl));
   free((profile  *) (m+nrl));
}




alignment   *avector(nl,nh)
 int nl, nh;
 /* Allocates an alignment vector with range [nl...nh].  */
 {
   alignment  *v;
   void nrerror();
   v=(alignment *)malloc((unsigned) (nh-nl+1)*sizeof(alignment));
   if (!v) nrerror("allocation failure in avector (check alignments dimension)");
   return v-nl;
 }

void free_avector(v,nl,nh)
alignment  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((alignment*) (v+nl));
}


sequence   *svector(nl,nh)
 int nl, nh;
 /* Allocates a sequence vector with range [nl...nh].  */
 {
   sequence  *v;
   void nrerror();
   v=(sequence *)malloc((unsigned) (nh-nl+1)*sizeof(sequence));
   if (!v) nrerror("allocation failure in svector");
   return v-nl;
 }

void free_svector(v,nl,nh)
sequence  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((sequence*) (v+nl));
}

reference   *rvector(nl,nh)
 int nl, nh;
 /* Allocates a sequence vector with range [nl...nh].  */
 {
   reference  *v;
   void nrerror();
   v=(reference *)malloc((unsigned) (nh-nl+1)*sizeof(reference));
   if (!v) nrerror("allocation failure in rvector");
   return v-nl;
 }

void free_rvector(v,nl,nh)
reference  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((reference*) (v+nl));
}


boundary   *bvector(nl,nh)
 int nl, nh;
 /* Allocates a sequence vector with range [nl...nh].  */
 {
   boundary  *v;
   void nrerror();
   v=(boundary *)malloc((unsigned) (nh-nl+1)*sizeof(boundary));
   if (!v) nrerror("allocation failure in bvector");
   return v-nl;
 }

void free_bvector(v,nl,nh)
boundary  *v;
int nl, nh;
/* Frees an sequence vector allocated by vector() */
{
 free((boundary*) (v+nl));
}


int   *ivector(nl,nh)
 int nl, nh;
 /* Allocates a sequence vector with range [nl...nh].  */
 {
   int  *v;
   void nrerror();
   v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
   if (!v) nrerror("allocation failure in ivector");
   return v-nl;
 }

void free_ivector(v,nl,nh)
int  *v;
int nl, nh;
/* Frees an sequence vector allocated by vector() */
{
 free((int*) (v+nl));
}


subalign  *sbvector(nl,nh)
 int nl, nh;
 /* Allocates an alignment vector with range [nl...nh].  */
 {
   subalign  *v;
   void nrerror();
   v=(subalign *)malloc((unsigned) (nh-nl+1)*sizeof(subalign));
   if (!v) nrerror("allocation failure in sbvector");
   return v-nl;
 }

void free_sbvector(v,nl,nh)
subalign  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((subalign*) (v+nl));
}

char  *cvector(nl,nh)
 int nl, nh;
 /* Allocates an alignment vector with range [nl...nh].  */
 {
   char  *v;
   void nrerror();
   v=(char *)malloc((unsigned) (nh-nl+1)*sizeof(char));
   if (!v) nrerror("allocation failure in cvector");
   return v-nl;
 }

void free_cvector(v,nl,nh)
char  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((char*) (v+nl));
}

float  *fvector(nl,nh)
 int nl, nh;
 /* Allocates an alignment vector with range [nl...nh].  */
 {
   float  *v;
   void nrerror();
   v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
   if (!v) nrerror("allocation failure in fvector");
   return v-nl;
 }

void free_fvector(v,nl,nh)
float  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((float*) (v+nl));
}

profile  *pvector(nl,nh)
 int nl, nh;
 /* Allocates a profile vector with range [nl...nh].  */
 {
   profile  *v;
   void nrerror();
   v=(profile *)malloc((unsigned) (nh-nl+1)*sizeof(profile));
   if (!v) nrerror("allocation failure in pvector");
   return v-nl;
 }

void free_pvector(p,nl,nh)
profile  *p;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((profile*) (p+nl));
}






