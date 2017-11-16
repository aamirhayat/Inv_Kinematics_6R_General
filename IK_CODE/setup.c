#include "function.h"
#include <stdio.h>
#include "global.h"

#define IMAGINARY 2.0e-3

int check_solutions(x3,x4,x5)
coordinate x3,x4,x5;

{
coordinate temp0[6][9];
coordinate temp1[6][9];
coordinate temp2[6][9];
coordinate temp3[6][9];
coordinate vector[9];
coordinate result;
int i,j;

vector[8] = 1;
vector[7] = x5;
vector[6] = x5 * x5;
vector[5] = x4* vector[8];
vector[4] = x4* vector[7];
vector[3] = x4* vector[6];
vector[2] = x4* vector[5];
vector[1] = x4* vector[4];
vector[0] = x4* vector[3];

for (i=0;i<6;i++)
 for (j=0;j<9;j++)
   temp3[i][j] = EM21[i][j] * x3*x3 + EM1[i][j]*x3 + EM0[i][j];

for (i=0;i<6;i++)
{ 
 result =0;
 for (j=0;j<9;j++)
   result += temp3[i][j] * vector[j];

   /* epslon is 0.1 here */
   if (fabs(result) > 0.22)
    return 0;
}

return 1;
}


 /* we know temp matrix are 6*6 matrix */

set_up_degen(sigma,verify_rank)
coordinate sigma[MAXG][MAXG];
int verify_rank;
{
 coordinate S_EM2[MAXD][MAXD],S_EM1[MAXD][MAXD],S_EM0[MAXD][MAXD];
 coordinate S_Dwr[MAXD2],S_Dwi[MAXD2],S_Dbeta[MAXD2],S_Dzz[MAXD2][MAXD2];
 coordinate temp0[MAXD][MAXD];
 coordinate temp1[MAXD][MAXD];
 coordinate temp2[MAXD][MAXD];
 coordinate temp3[MAXD][MAXD];
 coordinate temp[MAXD][MAXD];
 coordinate T_solution[24][3];
 int p[MAXG],q[MAXG];
 int i,j,k,k2,k3,rank;
 coordinate x3,x4,x5;

 double epslon;
 
 epslon = 1.0e-10;
 k3 =0;

 gauss(12,sigma,p,q);

 i = 11;

 while (fabs(sigma[p[i]][i]) < epslon)
   i -= 1;
 rank = i+1;
 
 for ( i=0;i< MAXD;i++)
 for ( j=0;j< MAXD;j++)
 {
  S_EM2[i][j] =0;
  S_EM1[i][j] =0;
  S_EM0[i][j] =0;
  }



 printf(" These are the ranks  RANK --%d , VERIFY --\n",rank);
 if (rank != verify_rank)
     printf(" THE RANK MAY NOT BE CORRECT!!!\n");
 else
     printf(" THE RANK IS CORRECT!!!\n");
    
 for (j=0;j<rank;j++)
  for (k=0;k<rank;k++)
   {
     S_EM2[j][k] = EM21[p[j]][q[k]];
     S_EM1[j][k] = EM1[p[j]][q[k]];
     S_EM0[j][k] = EM0[p[j]][q[k]];
    }
   reduce_eigen_degen(rank,S_EM2,S_EM1,S_EM0,S_Dwr,S_Dwi,S_Dbeta,S_Dzz);


   k2=0;
  for (k=0;k<2*rank;k++)
    if (fabs(S_Dwi[k])< IMAGINARY) 
     {
      D_wr[k2] = S_Dwr[k];
      D_wi[k2] = S_Dwi[k];
      k2+=1;
      }
  for (k=0;k<k2;k++)
    if (fabs(D_wi[k])< IMAGINARY) 
    {
    /*  if (!fGeneralized) */
 	for (i=0;i<6;i++)
	 for (j=0;j<3;j++)
	  {
	   S_EM2[i][j] = EM21[p[i]][j] * D_wr[k] * D_wr[k] + EM1[p[i]][j] * D_wr[k] + EM0[p[i]][j];
	   S_EM1[i][j] = EM21[p[i]][j+3] * D_wr[k] * D_wr[k] + EM1[p[i]][j+3] * D_wr[k] + EM0[p[i]][j+3];
	   S_EM0[i][j] = EM21[p[i]][j+6] * D_wr[k] * D_wr[k] + EM1[p[i]][j+6] * D_wr[k] + EM0[p[i]][j+6];
	   temp2[i][j] = S_EM2[i][j];
	   temp1[i][j] = S_EM1[i][j];
	   temp0[i][j] = S_EM0[i][j];
	   temp3[i][j] = temp2[i][j]+temp1[i][j]+temp0[i][j];
	   }
/*	   printf(" +++++++++++++++++++++++++++++++++++\n");
	   for (i=0;i<6;i++)
	    {
	     for (j=0;j<6;j++)
	      printf(" %f",temp2[i][j]);
	     printf("\n");
	      }
	   printf(" +++++++++++++++++++++++++++++++++++\n");
	   for (i=0;i<6;i++)
	    {
	     for (j=0;j<6;j++)
	      printf(" %f",temp1[i][j]);
	     printf("\n");
	      }
	   printf(" +++++++++++++++++++++++++++++++++++\n");
	   for (i=0;i<6;i++)
	    {
	     for (j=0;j<6;j++)
	      printf(" %f",temp0[i][j]);
	     printf("\n");
	      }
	   printf(" +++++++++++++++++++++++++++++++++++\n");
	   for (i=0;i<6;i++)
	    {
	     for (j=0;j<6;j++)
	      printf(" %f",temp3[i][j]);
	     printf("\n"); 
	      }
	      */

 for ( i=0;i< MAXD;i++)
 for ( j=0;j< MAXD;j++)
 {
  S_EM2[i][j] =0;
  S_EM1[i][j] =0;
  S_EM0[i][j] =0;
  }
	     for (j=0;j<3;j++)
	     {
	      S_EM2[0][j] = temp2[0][j];
	      S_EM2[1][j] = temp2[1][j];
	      S_EM1[0][j] = temp1[0][j];
	      S_EM1[1][j] = temp1[1][j];
	      S_EM0[0][j] = temp0[0][j];
	      S_EM0[1][j] = temp0[1][j];
	     }

	     for (j=0;j<3;j++)
	     {
	      S_EM2[2][j+1] = S_EM2[0][j];
	      S_EM2[3][j+1] = S_EM2[1][j];
	      S_EM1[2][j+1] = S_EM1[0][j];
	      S_EM1[3][j+1] = S_EM1[1][j];
	      S_EM0[2][j+1] = S_EM0[0][j];
	      S_EM0[3][j+1] = S_EM0[1][j];
	      }
	      


	   i = 4;
	reduce_eigen_degen(i,S_EM2,S_EM1,S_EM0,S_Dwr,S_Dwi,S_Dbeta,S_Dzz);
    
 
      printf("=========================\n");
      printf(" THis is the number of eigenvalue of X3[%d]  is %f\n",k,D_wr[k]);

      for (i=0;i<8;i++)
      {
    printf("value of X4 eigens--real %f    ---ima %f\n",S_Dwr[i],S_Dwi[i]);

      if (fabs(S_Dwi[i]) < IMAGINARY)
       {
	x3= D_wr[k];
	x4= S_Dwr[i];
	x5 = S_Dzz[i][0]/S_Dzz[i][1];
	
	printf(" x3 = %f, x4= %f, x5 =%f\n",x3,x4,x5);
       if (check_solutions(x3,x4,x5))
	{
	 T_solution[k3][0] = x3;
	 T_solution[k3][1] = x4;
	 T_solution[k3][2] = x5;
	 k3 += 1;
	 }
       }
      }
   }

   for (i=0;i<k3;i++)
    printf(" Temp solutions x3 %f x4 %f x5 %f\n",T_solution[i][0],T_solution[i][1],T_solution[i][2]);

  for (i=0;i<24;i++)
  {
   D_wr[i]=0;
   D_wi[i] = 1;
   D_beta[i] =1;
   }
   for (i=0;i<k3;i++)
    {
     D_wr[i] = T_solution[i][0];
     D_wi[i] =0;
     D_beta[i] =1;
     fGeneralized=0;
     D_zz[i][0] = T_solution[i][1] * T_solution[i][2];
     D_zz[i][1] = T_solution[i][1];
     D_zz[i][3] = T_solution[i][2];
     }
}
