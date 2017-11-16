#include <stdio.h>
#include "function.h"
#include "global.h"

int check_degen(sigma,tag)
coordinate sigma[MAXG][MAXG];
int *tag;
{
 coordinate total[MAXG][MAXG];
 coordinate r_x3;
 int p[MAXG],q[MAXG];
 //FILE *fp1,*fp2,*fp3,*fp4,*fp5;

 
 int rank,rank2;
 int i,j,k;
/*
 fp1 = fopen("check1.dat","w");
 fp2 = fopen("check2.dat","w");
 fp3 = fopen("sigma2.dat","w");
 fp4 = fopen("sigma1.dat","w");
 fp5 = fopen("sigma0.dat","w");
*/
 r_x3 = 0.4;

    for (i=0;i<MAXG;i++)
     for (j=0;j<MAXG;j++)
     {
      total[i][j] = 0;
     sigma[i][j] =0;
     }

  for (i=0;i<12;i++)
   for (j=0;j<12;j++)
     sigma[i][j] = EM21[i][j] * r_x3 * r_x3 + EM1[i][j] *r_x3 + EM0[i][j];


/*   for (i=0;i<9;i++)
   {
    for (j=0;j<9;j++)
     {
       fprintf(fp1,"%f ",sigma[i][j]);
       fprintf(fp3,"%f ",EM21[i][j]);
       fprintf(fp4,"%f ",EM1[i][j]);
       fprintf(fp5,"%f ",EM0[i][j]);
       }
    fprintf(fp1,"\n");
    fprintf(fp3,"\n");
    fprintf(fp4,"\n");
    fprintf(fp5,"\n");
   } */


    for (i=0;i<9;i++)
    {
     for (j=0;j<6;j++)
      total[i][j] = EM21[j][i];
     for (j=6;j<12;j++)
      total[i][j] = EM1[j-6][i];
     for (j=12;j<18;j++)
      total[i][j] = EM0[j-12][i];
    }
  /* for (i=0;i<18;i++)
   {
    for (j=0;j<18;j++)
       fprintf(fp2,"%f ",total[i][j]);
    fprintf(fp2,"\n");
   }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5); */

  gauss(12,sigma,p,q);

  i = 11;
  while (fabs(sigma[p[i]][i]) < 0.00001)
   i -= 1;
  
  rank = i+1;

  for (i=0;i<12;i++)
   for (j=0;j<12;j++)
     sigma[i][j] = EM21[i][j] * r_x3 * r_x3 + EM1[i][j] *r_x3 + EM0[i][j];


  if (rank == 12)
  {
  *tag =0;
  return 0;
  }
   
  printf(" This is the first rank %d\n",rank);
  gauss(18,total,p,q);

  i = 17;
  while (fabs(total[p[i]][i]) < 0.00001)
   i -= 1;
   rank2 = i+1;

  printf(" This is the Second rank %d\n",rank2);

  if (rank == 2*rank2)
   *tag =1;
  else 
   *tag =2;

   return rank;


}


