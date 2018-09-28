#include  "inc.h"

void PrintPSSM(len_a,len_b,matrix)
int   len_b, len_a;
float **matrix;
{
  int   i,j,k;

  printf("*** PRINT DATA ***:\n");
  printf("        ");
  for(j=1;j<=len_b;j++)
  {
    printf("%8d|",j);
  }
  printf("\n");
  for(j=1;j<=len_a;j++)
  {
   printf("--------|");
   for(i=1;i<=len_b;i++){printf("---------");}
   printf("\n");
   printf("%8d|",j);
   for(i=1;i<=len_b;i++)
   {
    printf("%8.1e;",matrix[i][j]);
   }
   printf("\n");
  }
  printf("---------|");
  for(i=1;i<=len_b;i++)
  {
   printf("----------+");
  }
  printf("\n");


}
