#include  "inc.h"

void PrintMatrix(len_a,len_b,dim_data,data_matrix,data_seq)
int   len_b, len_a, dim_data;
float data_matrix[10][MAXS][MAXS];
int   data_seq[2][MAXS];
{
  int   i,j,k;

  printf("*** PRINT DATA ***:\n");
  printf("          ");
  for(j=0;j<len_b;j++)
  {
    for (k=0;k<dim_data-1;k++){printf("         ");}
    printf("  %8d|",data_seq[1][j]);
  }
  printf("\n");
  for(j=0;j<len_a;j++)
  {
   printf("---------|");
   for(i=0;i<len_b;i++){for (k=0;k<dim_data-1;k++){printf("----------");}}
   printf("\n");
   printf(" %8d|",data_seq[0][j]);
   for(i=0;i<len_b;i++)
   {
    printf("(");
    for (k=0;k<dim_data-1;k++){
     if (data_matrix[k][i][j] >=0.0) {printf("%8.1e;",data_matrix[k][i][j]);
     }else{printf("%8.1e;",data_matrix[k][i][j]);}
    }
    if (data_matrix[dim_data-1][i][j] >=0.0) {printf("%8.1e)|",data_matrix[dim_data-1][i][j]);
    }else{printf("%8.1e;",data_matrix[dim_data-1][i][j]);}
  }
   printf("\n");
  }
  for(i=0;i<len_b;i++)
  {
   printf("----------");
  }
  printf("----------\n");
  printf("\n");


}
