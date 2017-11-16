#include "function.h"
#include "global.h"
#include <stdio.h>
#include <stdlib.h>

/* 
  Given EM2, EM1, EM0, it reduces the univariate equation solving to an 
  eigenvalue problem.
*/

reduce_Eigen()
{
   int i,j,k,i1;
   int dim,order, ipvt[12], job;
   coordinate z[12];
   coordinate det1[2];
   coordinate fv1[24];
   int ierr, matz, iv1[24];
   coordinate EPS = 0.00001;
   coordinate temp1,temp2,sum, sum1;
   coordinate rshift,cshift;
   extern double sqrt();
   coordinate eps1 = 0.0000005;
   char balanc,compvl,compvr,sense;
   int ihi,ilo,info,lda,ldvl,ldvr,lwork,iwork[200];
   char compq;
   int select[24],M,index[24];
   coordinate S,SEP,WR[24],WI[24],WORK[200],SUM,U[12][12], V[12][12], SS[12], E[12];
   int degenerated;
   int verify_rank;

   
   verify_rank = 12;
   degenerated = 0;
   matz = 10;
   fGeneralized = 0;
   dim = 12;
   for (i = 0; i < 12; i++)
   {
       ipvt[i] = 0;
       D_beta[i] = D_beta[i+12] = 1.0;
   }

   Icond = 0.0;
   det1[0] = det1[1] = 0.0;
   dgeco_(EM2,&dim,&dim,ipvt,&Icond,z);
   printf(" Icond = %f \n",Icond);
   if (Icond >= EPS)
   {
    
       job = 11;
       dgedi_(EM2,&dim,&dim,ipvt,det1,work,&job);
       Mat_Mul(EM2,EM1,IEM1,dim);
       Mat_Mul(EM2,EM0,IEM0,dim);
    
       for (i = 0; i < 12; i++)
       for (j = 0; j < 24; j++)
           EIG[j][i] = 0.0;
       
       for (i = 0; i < 12; i++)
           EIG[i+12][i] = 1.0;
    
       for (i = 12; i < 24; i++)
       for (j = 0; j < 12; j++)
           EIG[j][i] = - IEM0[i-12][j];
    
       for (i = 12; i < 24; i++)
       for (j = 12; j < 24; j++)
           EIG[j][i] = - IEM1[i-12][j-12];
    
       dim = 24;
       
       for (i = 0; i < 24; i++)
       for (j = 0; j < 24; j++)
	EIGcopy[i][j]= EIG[i][j];

       rg_(&dim,&dim,EIG,D_wr,D_wi,&matz,D_zz,iv1,fv1,&ierr);

       for (i = 0; i < 24; i++)
           printf("EIG: %d %lf %lf \n",i, D_wr[i],D_wi[i]);


       if (ierr)
          printf(" Eigenvalue Computation : ierr = %f \n", ierr);
             
    }
    else
    {  
        if (Perform_Trans())
        {
            Mat_Mul(NEM2,NEM1,NIEM1,dim);
            Mat_Mul(NEM2,NEM0,NIEM0,dim);
        
            for (i = 0; i < 12; i++)
            for (j = 0; j < 24; j++)
                EIG[j][i] = 0.0;
            
            for (i = 0; i < 12; i++)
                EIG[i+12][i] = 1.0;
        
            for (i = 12; i < 24; i++)
            for (j = 0; j < 12; j++)
                EIG[j][i] = - NIEM0[i-12][j];
            
            for (i = 12; i < 24; i++)
            for (j = 12; j < 24; j++)
                EIG[j][i] = - NIEM1[i-12][j-12];
        
            dim = 24;

            rg_(&dim,&dim,EIG,D_wr,D_wi,&matz,D_zz,iv1,fv1,&ierr);

            for (i = 0; i < 24; i++)
            {
               temp1 = D_wr[i];
               temp2 = D_wi[i];
/*
   Now check whether the original system had an eigenvalue close to infinity.
   In the transformed system that would correspond to an eigenvalue very
   close to -Rdelta/Rgamma.
*/
               sum = Rgamma * temp1 + Rdelta; 
               sum = sum * sum;
               sum = sum + Rgamma * temp2 * Rgamma * temp2;
               sum1 = sqrt(sum);
               if (sum < eps1)
               {
                       D_wi[i] = 0.0;
                       D_wr[i] = 1.0e64;
               }
               else
               {
               D_wr[i] = ((Ralpha * temp1 + Rbeta) * (Rgamma * temp1 + Rdelta) + Ralpha * Rgamma * temp2 * temp2) / sum;
               D_wi[i] =  temp2 * (Ralpha * Rdelta - Rgamma * Rbeta) / sum;
               }
               printf("LINEAR TRANSFORM: %d. %15.10g + i %15.10g, %15.10g + i %15.10g \n",i,temp1, temp2, D_wr[i], D_wi[i]);
           }
        }

     /* We solve it as a Generalized Eigenvalue Problem 
               of the form EIG1 * q + EIG   
      */
  else
        {
           fGeneralized = 1;
           for (i = 0; i < 12; i++)
           for (j = 0; j < 24; j++)
               TT[j][i] = EIG[j][i] = 0.0;
           
           for (i = 0; i < 12; i++)
               TT[i+12][i] = EIG[i+12][i] = 1.0;
        
           for (i = 12; i < 24; i++)
           for (j = 0; j < 12; j++)
               TT[j][i] = EIG[j][i] = -EM0[i-12][j];
    
           for (i = 12; i < 24; i++)
           for (j = 12; j < 24; j++)
               TT[j][i] = EIG[j][i] = -EM1[i-12][j-12];
        
           for (i = 0; i < 12; i++)
           for (j = 0; j < 24; j++)
               TT1[j][i] = EIG1[j][i] = 0.0;
           
           for (i = 0; i < 12; i++)
               TT1[i][i] = EIG1[i][i] = 1.0;
    
           for (i = 12; i < 24; i++)
           for (j = 0; j < 12; j++)
               TT1[j][i] = EIG1[j][i] = 0.0;
        
           for (i = 12; i < 24; i++)
           for (j = 12; j < 24; j++)
               TT1[j][i] = EIG1[j][i] = EM21[i-12][j-12];

           for (i = 0; i < 24; i++)
           for (j = 0; j < 24; j++)
	    {
	     EIGcopy[i][j] = EIG[i][j];
	     EIG1copy[i][j] = EIG1[i][j];
	     }

           
           dim = 24;
    /*
    
          The matrix equation system has an eigenvalue at infinity, iff the leading
          matrix, EG21, is singular. We may check that separately.
    
    */

           rgg_(&dim,&dim,EIG,EIG1,D_wr,D_wi,D_beta,&matz,D_zz,&ierr);

           for (i = 0; i < 24; i++)
               printf("general EIG: %d %lf %lf \n",i, D_wr[i]/D_beta[i],D_wi[i]/D_beta[i]);

           for (i = 0; i < 24; i++)
               if (D_beta[i] != 0.0)
               {
		  verify_rank = 12;
                   D_wr[i] = D_wr[i]/D_beta[i] ;
                   D_wi[i] = D_wi[i]/D_beta[i] ;
                   D_beta[i] = 1.0 ;
                   if (fabs(D_wi[i]) <= 1.0e-2/Rcond)
                       printf("D_wr[i] = %17.13f \n",D_wr[i]);
		}
	}


/*		for (i=0;i<24;i++)
                   if (fabs(D_wi[i]) <= 1.0e-2/Rcond)
		   {
                       for(j = 0; j < 12; j++)
                          for(k = 0; k < 12; k++)
                              SVD1[k][j] = SVD[j][k] = EM21[j][k] * D_wr[i] * D_wr[i] + EM1[j][k] * D_wr[i] + EM0[j][k];


                      dim = 12;
                      job = 01;
                      dsvdc_(SVD1,&dim,&dim,&dim,SS,E,U,&dim,V,&dim,WORK,&job,&ierr);
               if (SS[11] > 1.0e-6)
                  printf("gen Eigenvalue is not Accurate \n");
               if (SS[10] < 1.0e-11)
		 { degenerated = 1;
                  printf("gen There are Possibly INFINITE SOLUTIONS \n");
		  }
		  while (SS[verify_rank-1] < 1.0e-9)
		   verify_rank -= 1;

                      for (j = 0; j < 12; j++)
                      {
                          D_zz[i][j] = V[11][j];
                          D_zz[i][j+12] = V[11][j] * D_wr[i];
                      }
		     }
*/

   }
}

