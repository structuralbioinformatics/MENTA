#include        <stdio.h>
#include        <math.h>
#include        <ctype.h>
#include        <string.h>
#include        <sys/types.h>
#include        <malloc.h>
#define         MAXS  500   /* minimum 300 */
#define         MAXG 100   
#define         DIM  1000
#define         MAXV  15.0 
#define         MAXVG 100.0 
#define         MAXL  100.0
#define         MAXQ  700
#define         MAXF  100
#define         MAXFS 700
#define         MAXP  10
#define         MAXM  10
#define         MAXE  30
#define         MAXSC 1000.0
#define         MAXTI 15
#define         MAXTITLE 2000
#define         MAXLINE 80
#define         MAXMETH 14
#define         MAXEVD 500
#define         MIN   -1.0e+12
#define         MINP   1.0e-5
#define         DATASIZE 5000
#define         PMURZIN  0.1
#define         MAXTOL   20
#define         MAXPROF  1000
#define         MAXPS    100

int EVAL;

typedef struct consensus
{
 char     feature[5];
} consensus;

typedef struct sequence
{
 char    title[MAXTITLE];        /* INPUT TITLES HAVE LESS THAN 50 CHARACTERS */
 int     element[MAXS];
 int     position[MAXS];
 int     length;
} sequence;

typedef struct reference
{
 int     image[MAXS];
 int     position[MAXS];
} reference;

typedef struct alignment
{
 char      title[MAXTITLE];    /*  OUTPUT TITLES HAVE LESS THAN 1000  CHARACTERS */
 sequence  template;
 sequence  target;
 consensus profile[MAXS];
 reference Itemplate;
 reference Itarget;
 int       length;
 float     score;
 double    evalue;
 double    pvalue;
 float     homol;
 float     ident;
} alignment;


typedef struct profile
{
 int          cluster_dimension;
 int          cluster[MAXPROF];
 int          size;
 int          align;
 int          action;
 float        score;
 double       evalue;
 double       pvalue;
 float        homol;
 float        ident;
 sequence     main;
 sequence     sequence[MAXFS];
 consensus    profile[MAXS];
} profile;

typedef struct tabulation
{
  float   value;
  float   score;
  int     up;   
  int     left; 
  int     match;  
  int     backu;   
  int     backl;    
  int     dimension;      
} tabulation;

typedef struct subalign
{
  int     group;
  int     Ipair[MAXG];
  int     Jpair[MAXG];
  int     Imax;
  int     Jmax;
} subalign;


typedef struct boundary
{
  int     imax;
  int     jmax;
  int     imin;
  int     jmin;
} boundary;


typedef struct max3
{
  float   value;
  float   score;
  int     first;
  int     second;
  int     third;
} max3;

typedef struct psic
{
  float   element[MAXE];
  float   frequency[MAXE];
  float   targetQ[MAXE];
} psic;

typedef struct evd
{
  float   K;
  float   lambda;
  float   sigma;
  int     Rc;
  float   Sc;
  float   min;
} evd;

/*  PROGRAM  MENTA:  Multiple Entities Alignment */
/*  Baldomero Oliva. Computational Structural Biology Laboratory */
/*  Universitat Pompeu Fabra */
/*  Barcelona. Catalonia, Spain (EU */

