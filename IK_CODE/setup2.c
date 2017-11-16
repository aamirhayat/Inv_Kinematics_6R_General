#include "function.h"
#include <stdio.h>
#include "global.h"

#define IMAGINARY 2.0e-3

set_up_degen_32(sigma,verify_rank)
coordinate sigma[MAXG][MAXG];
int verify_rank;
{
 coordinate S_EM2[MAXD][MAXD],S_EM1[MAXD][MAXD],S_EM0[MAXD][MAXD];
 coordinate S_Dwr[MAXD2],S_Dwi[MAXD2],S_Dbeta[MAXD2],S_Dzz[MAXD2][MAXD2];
 coordinate temp0[12][13];
 coordinate temp1[12][13];
 coordinate temp2[12][13];
 coordinate temp3[12][12];
 coordinate temp[12][12];
 int p[MAXG],q[MAXG],s[12];
 int i,j,k,rank,count;

 double epslon;
 
 epslon = 1.0e-10;
/* rank = 2* rank; */
 gauss(12,sigma,p,q);
 i = 11;

 while (fabs(sigma[p[i]][i]) < epslon)
   i -= 1;
 rank = i+1;
 

 for ( i=0;i< MAXD;i++)
 for ( j=0;j<MAXD ;j++)
 {
  S_EM2[i][j] =0;
  S_EM1[i][j] =0;
  S_EM0[i][j] =0;
  }

#ifndef IK_OPT

/*
 printf(" These are the ranks  RANK --%d , VERIFY --\n",rank);

 if (rank != verify_rank)
     printf(" THE RANK MAY NOT BE CORRECT!!!\n");
 else
     printf(" THE RANK IS CORRECT!!!\n");


 for (j=0;j<10;j++)
 printf(" This is p[%d] %d\n",j,p[j]); */

#endif


 for (i=0;i<12;i++)
 for(j=0;j<13;j++)
  {
   temp2[i][j] =0;
   temp1[i][j] =0;
   temp0[i][j] =0;
   }

    
j=0;

  for (count=0;count<4;count++)
  {
      while (p[j]>=6) j++;
      if (j >6 )
	    printf(" THERE ARE POSSIBLY INFINITE SOLUTIONs, quit\n");

  for (k=0;k<9;k++)
   {
     temp2[count][k] = EM21[p[j]][k];
     temp1[count][k] = EM1[p[j]][k];
     temp0[count][k] = EM0[p[j]][k];
    }
    j++;
  }

 for (j=0;j<4;j++)
  for (k=0;k<9;k++)
   {
    temp2[j+4][k+3] = temp2[j][k];
    temp1[j+4][k+3] = temp1[j][k];
    temp0[j+4][k+3] = temp0[j][k];
    }


 s[0]=0; s[1]=3; s[2]=6; s[3]=9; s[4]=1;
 s[5]=2; s[6]=12; s[7]=4; s[8]=5; s[9]=12;
 s[10]=7; s[11]=8; s[12]=12; s[13]=10; s[14]=11; s[15]=12;


 for (i=0;i<8;i++)
  for (j=0;j<16;j++)
   {
  S_EM2[i][j] = temp2[i][s[j]];
  S_EM1[i][j] = temp1[i][s[j]];
  S_EM0[i][j] = temp0[i][s[j]];
  }

  for (i=0;i<8;i++)
   for (j=0;j<12;j++)
   {
  S_EM2[i+8][j+4] = temp2[i][j];
  S_EM1[i+8][j+4] = temp1[i][j];
  S_EM0[i+8][j+4] = temp0[i][j];
  }
    
rank = 16;
    
   reduce_eigen_degen(rank,S_EM2,S_EM1,S_EM0,S_Dwr,S_Dwi,S_Dbeta,S_Dzz);

   for (i=0;i<24;i++)
   {
    D_wr[i] = 0;
    D_wi[i] = 1;
    D_beta[i] = 1;
    fGeneralized =0;
    }

j=0;
for (i=0;i<MAXD2;i++)
 {

 /* 0.2 is a randomly chosed esplon */
  if (fabs(S_Dwi[i])< 0.2)
  {
   D_wr[j] = S_Dwr[i];
   D_wi[j] = S_Dwi[i];
   D_beta[j] = S_Dbeta[i];

   for (k=0;k<12;k++)
    D_zz[j][k] = S_Dzz[i][k+4];
   for (k=12;k<24;k++)
    D_zz[j][k] = S_Dzz[i][k+8];
    j+=1;
   }
  }
  
  if (j>12) printf(" To many real solutions\n");

 
}

