#include "function.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "global.h"

#define MAX_TRANS 3  /* Perform_trans_degen use it constant */

/* 
  Given EM2, EM1, EM0, it reduces the univariate equation solving to an 
  eigenvalue problem.
*/

reduce_eigen_degen(num,EM2,EM1,EM0,D_wr,D_wi,D_beta,D_zz)
int num;
coordinate EM2[MAXD][MAXD],EM1[MAXD][MAXD],EM0[MAXD][MAXD];
coordinate D_wr[32], D_wi[32],D_beta[32],D_zz[32][32];

{
  coordinate NEM2[MAXD][MAXD],NEM21[MAXD][MAXD];
  coordinate NEM1[MAXD][MAXD],INEM1[MAXD][MAXD];
  coordinate NEM0[MAXD][MAXD],INEM0[MAXD][MAXD]; 

  coordinate IEM1[MAXD][MAXD],IEM0[MAXD][MAXD]; 


  coordinate EIG[MAXD2][MAXD2], EIG1[MAXD2][MAXD2],EIGcopy[MAXD2][MAXD2],EIG1copy[MAXD2][MAXD2];
  coordinate TT[MAXD2][MAXD2], TT1[MAXD2][MAXD2];
  coordinate SVD1[MAXD][MAXD],SVD[MAXD][MAXD];  

   FILE *fp;
   int i,j,k,i1;
   int dim,order, ipvt[MAXD], job;
   coordinate z[MAXD];
   coordinate det1[2];
   coordinate fv1[MAXD2];
   int ierr, matz, iv1[MAXD2];
   coordinate EPS = 0.0000001;
   coordinate temp1,temp2,sum, sum1;
   coordinate rshift,cshift;
   extern double sqrt();
   coordinate eps1 = 0.0000005;
   char balanc,compvl,compvr,sense;
   int ihi,ilo,info,lda,ldvl,ldvr,lwork,iwork[200];
   char compq;
   int select[MAXD2],M,index[MAXD2];
   coordinate S,SEP,WR[MAXD2],WI[MAXD2],WORK[200],SUM,U[MAXD][MAXD], V[MAXD][MAXD], SS[MAXD], E[MAXD2];
   coordinate EM21[MAXD][MAXD];
   coordinate ver_vector[MAXD2];
   coordinate MAT_NUM[1024],MAT_NUM1[1024];
   coordinate MAT_NUM2[1024],MAT_NUM21[1024],MAT_NUM22[1024];
   int num2;
   int it1,jt1;


   
   num2= 2*num;

   
   for (i = 0; i < num; i++)
     for (j = 0; j < num; j++)
       EM21[i][j] = EM2[i][j];

   matz = 10;
   fGeneralized = 0;
   dim = num;

   for (i = 0; i < num; i++)
   {
       ipvt[i] = 0;
       D_beta[i] = D_beta[i+num] = 1.0;
   }

   Icond = 0.0;

   det1[0] = det1[1] = 0.0;

   for (i = 0; i < num; i++)
     for (j = 0; j < num; j++)
       MAT_NUM[i*num+j] = EM2[j][i];  
   dgeco_(MAT_NUM,&dim,&dim,ipvt,&Icond,z);
   for (i = 0; i < num; i++)
     for (j = 0; j < num; j++)
       EM2[j][i] = MAT_NUM[i*num+j];

   printf(" Icond = %f \n",Icond);

   if (Icond >= EPS)
   {
    
       job = 11;
       for (i = 0; i < num; i++)
         for (j = 0; j < num; j++)
           MAT_NUM[i*num+j] = EM2[j][i];  
       dgedi_(MAT_NUM,&dim,&dim,ipvt,det1,work,&job);
       for (i = 0; i < num; i++)
         for (j = 0; j < num; j++)
           EM2[j][i] = MAT_NUM[i*num+j];

       Mat_MulD(EM2,EM1,IEM1,dim);
       Mat_MulD(EM2,EM0,IEM0,dim);
    
       for (i = 0; i < num; i++)
       for (j = 0; j < num2; j++)
           EIG[j][i] = 0.0;
       
       for (i = 0; i < num; i++)
           EIG[i+num][i] = 1.0;
    
       for (i = num; i < num2; i++)
       for (j = 0; j < num; j++)
           EIG[j][i] = - IEM0[i-num][j];
    
       for (i = num; i < num2; i++)
       for (j = num; j < num2; j++)
           EIG[j][i] = - IEM1[i-num][j-num];
    
       dim = num2;
       
       for (i = 0; i < num2; i++)
       for (j = 0; j < num2; j++)
	EIGcopy[i][j]= EIG[i][j];

       for (i = 0; i < num2; i++)
	  for (j = 0; j < num2; j++)
	    {
	     MAT_NUM2[i*num2+j] = EIG[i][j];
	     MAT_NUM21[i*num2+j] = D_zz[j][i];
	     }
       rg_(&dim,&dim,MAT_NUM2,D_wr,D_wi,&matz,MAT_NUM21,iv1,fv1,&ierr);
       for (i = 0; i < num2; i++)
	  for (j = 0; j < num2; j++)
	   {
	     EIG[i][j] = MAT_NUM2[i*num2+j];
	     /* order reversing has been considered */
	     D_zz[i][j] = MAT_NUM21[i*num2+j] ;
	     }
/*

       for (i = 0; i < num2; i++)
       {
	  printf("EIG[%d] = %.12lf %.12lf \n",i, D_wr[i],D_wi[i])
	  if (num2 != 32)
	  {
	  printf(" EIG for  x4 is %.12lf\n", D_zz[i][0]/D_zz[i][3]); 
	  printf(" EIG for  x5 is %.12lf\n", D_zz[i][0]/D_zz[i][1]); 
	  }
	  if (num2 == 32)
	  {
	  printf(" EIG for  x4 is %.12lf\n", D_zz[i][4]/D_zz[i][7]); 
	  printf(" EIG for  x5 is %.12lf\n", D_zz[i][4]/D_zz[i][5]); 
	  }
	  printf("\n");
	 }


*/
       if (ierr)
          printf(" Eigenvalue Computation : ierr = %f \n", ierr);
             
   }
    else
    {  
    /* the Perform_Trans_degen has not been test yet */
      if (Perform_Trans_degen(num,EM21,EM1,EM0,NEM2,NEM1,NEM0))
        {
            Mat_MulD(NEM2,NEM1,NIEM1,dim);
            Mat_MulD(NEM2,NEM0,NIEM0,dim);
        
            for (i = 0; i < num; i++)
            for (j = 0; j < num2; j++)
                EIG[j][i] = 0.0;
            
            for (i = 0; i < num; i++)
                EIG[i+num][i] = 1.0;
        
            for (i = num; i < num2; i++)
            for (j = 0; j < num; j++)
                EIG[j][i] = - NIEM0[i-num][j];
            
            for (i = num; i < num2; i++)
            for (j = num; j < num2; j++)
                EIG[j][i] = - NIEM1[i-num][j-num];
        
            dim = num2;
            for (i = 0; i < num2; i++)
              for (j = 0; j < num2; j++)
               {
                MAT_NUM2[i*num2+j] = EIG[i][j];
                MAT_NUM21[i*num2+j] = D_zz[i][j];
                }
            rg_(&dim,&dim,MAT_NUM2,D_wr,D_wi,&matz,MAT_NUM21,iv1,fv1,&ierr);
            for (i = 0; i < num2; i++)
              for (j = 0; j < num2; j++)
               {
                EIG[i][j] = MAT_NUM2[i*num2+j];
                D_zz[i][j] = MAT_NUM21[i*num2+j];
                }
            for (i = 0; i < num2; i++)
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
            /*   printf("LINEAR TRANSFORM: %d. %15.10g + i %15.10g, %15.10g + i %15.10g \n",i,temp1, temp2, D_wr[i], D_wi[i]);
	       printf(" --------This is is the eigenvectors \n");
	       for (k=0; k< num2;k++)
	 	printf(" D_zz[%d] == %f\n",k,D_zz[i][k]);  */
           }
        }

     /* We solve it as a Generalized Eigenvalue Problem 
               of the form EIG1 * q + EIG   
      */
      else
        {
           fGeneralized = 1;
           for (i = 0; i < num; i++)
           for (j = 0; j < num2; j++)
               TT[j][i] = EIG[j][i] = 0.0;
           
           for (i = 0; i < num; i++)
               TT[i+num][i] = EIG[i+num][i] = 1.0;
        
           for (i = num; i < num2; i++)
           for (j = 0; j < num; j++)
               TT[j][i] = EIG[j][i] = -EM0[i-num][j];
    
           for (i = num; i < num2; i++)
           for (j = num; j < num2; j++)
               TT[j][i] = EIG[j][i] = -EM1[i-num][j-num];
        
           for (i = 0; i < num; i++)
           for (j = 0; j < num2; j++)
               TT1[j][i] = EIG1[j][i] = 0.0;
           
           for (i = 0; i < num; i++)
               TT1[i][i] = EIG1[i][i] = 1.0;
    
           for (i = num; i < num2; i++)
           for (j = 0; j < num; j++)
               TT1[j][i] = EIG1[j][i] = 0.0;
        
           for (i = num; i < num2; i++)
           for (j = num; j < num2; j++)
               TT1[j][i] = EIG1[j][i] = EM21[i-num][j-num];

           for (i = 0; i < num2; i++)
           for (j = 0; j < num2; j++)
	    {
	     EIGcopy[i][j] = EIG[i][j];
	     EIG1copy[i][j] = EIG1[i][j];
	     }

           
           dim = num2;
    /*
    
          The matrix equation system has an eigenvalue at infinity, iff the leading
          matrix, EG21, is singular. We may check that separately.
    
    */
	   for (i = 0; i < num2; i++)
	      for (j = 0; j < num2; j++)
	       {
		 MAT_NUM2[i*num2+j] = EIG[j][i];
		 MAT_NUM21[i*num2+j] = EIG1[j][i];
		 MAT_NUM22[i*num2+j] = D_zz[j][i];
		 }
           rgg_(&dim,&dim,MAT_NUM2,MAT_NUM21,D_wr,D_wi,D_beta,&matz,MAT_NUM22,&ierr);
	   for (i = 0; i < num2; i++)
	      for (j = 0; j < num2; j++)
	      {
		 EIG[j][i]= MAT_NUM2[i*num2+j];
		 EIG1[j][i]= MAT_NUM21[i*num2+j];
		 D_zz[j][i]= MAT_NUM22[i*num2+j];
		 }
           for (i = 0; i < num2; i++)
              /* printf("general EIG: %d %lf %lf \n",i, D_wr[i]/D_beta[i],D_wi[i]/D_beta[i]);
*/
           for (i = 0; i < num2; i++)
	      if (D_beta[i] != 0.0)
		 {
		     D_wr[i] = D_wr[i]/D_beta[i] ;
		     D_wi[i] = D_wi[i]/D_beta[i] ;
		     D_beta[i] = 1.0 ;
		   /*  if (fabs(D_wi[i]) <= 1.0e-2/Rcond)
			  printf("D_wr[i] = %17.13f \n",D_wr[i]); */
		  }
	}	
           for (i=0;i<10;i++)
           {
              
               for(j = 0; j < num; j++)
                  {
                      for(k = 0; k < num; k++)
                         {

                             SVD1[k][j] = SVD[j][k] = EM21[j][k] * i*i + EM1[j]
[k] * i + EM0[j][k];
                         }
                  }
                      dim = num;
                      job = 01;

                          for (it1 = 0; it1 < num; it1++)
                             for (jt1 = 0; jt1 < num; jt1++)
                              {
                                MAT_NUM[it1*num+jt1]= SVD1[jt1][it1];
                                MAT_NUM1[it1*num+jt1]= V[jt1][it1];
                                }

               dsvdc_(MAT_NUM,&dim,&dim,&dim,SS,E,U,&dim,MAT_NUM1,&dim,WORK,&job,&ierr);
                          for (it1 = 0; it1 < num; it1++)
                             for (jt1 = 0; jt1 < num; jt1++)
                             {
                                SVD1[jt1][it1] = MAT_NUM[it1*num+jt1];
                                V[jt1][it1] = MAT_NUM1[it1*num+jt1];
                                }


              /* printf(" For X3 = %d \n", i);

               for (j=0;j<num;j++)
                         printf("gen SS[ %d ] --> %.15lf\n",j,SS[j]); */
           }




		for (i=0;i<num2;i++)
		     if (fabs(D_wi[i]) <= 1.0e-2/Rcond)
		      {
			  for(j = 0; j < num; j++)
			     for(k = 0; k < num; k++)
				SVD1[k][j] = SVD[j][k] = EM21[j][k] * D_wr[i] * D_wr[i] + EM1[j][k] * D_wr[i] + EM0[j][k];


			  dim = num;
			  job = 01;
			  for (it1 = 0; it1 < num; it1++)
			     for (jt1 = 0; jt1 < num; jt1++)
			     {
				MAT_NUM[it1*num+jt1] = SVD1[jt1][it1];  
				MAT_NUM1[it1*num+jt1] = V[jt1][it1];  
				}
			  dsvdc_(MAT_NUM,&dim,&dim,&dim,SS,E,U,&dim,MAT_NUM1,&dim,WORK,&job,&ierr);
			  for (it1 = 0; it1 < num; it1++)
			     for (jt1 = 0; jt1 < num; jt1++)
			     {
				SVD1[jt1][it1] = MAT_NUM[it1*num+jt1];  
				V[jt1][it1] = MAT_NUM1[it1*num+jt1];  
				}
			  if (SS[num-1] > 1.0e-6)
			     printf("gen Eigenvalue is not Accurate \n");
			  if (SS[num-2] < 1.0e-11)
			     printf("gen There are Possibly INFINITE SOLUTIONS \n");

			  for (j = 0; j < num; j++)
			     {
				 D_zz[i][j] = V[num-1][j];
				 D_zz[i][j+num] = V[num-1][j] * D_wr[i];
			     }

		      }
	}
}



int Perform_Trans_degen(num,EM21,EM1,EM0,NEM2,NEM1,NEM0)
int num;
coordinate EM21[MAXD][MAXD],EM1[MAXD][MAXD],EM0[MAXD][MAXD];
coordinate NEM2[MAXD][MAXD],NEM1[MAXD][MAXD],NEM0[MAXD][MAXD];

{
    coordinate EPS ;
    int i,j,k;
    coordinate temp1, temp2;
    int seed;
    coordinate random, z[MAXD];
    extern float urand_();
    int count, dim, ipvt[MAXD], job;
    coordinate det1[2];
    coordinate work[MAXD];
    coordinate factor;
    coordinate MAT_NUM[1024];
    EPS = 0.00002;


    seed = 10;
    count = 0;
    dim = num;
    job = 11;
    factor = 1.0;
    srand(seed);
    while ( count < MAX_TRANS)
    { 
        random = ( (double) rand() / RAND_MAX) ;
        if (random < 0.5)
            random = -1.0 + random;
        Ralpha = random;

        random = ( (double) rand() / RAND_MAX) ;
        if (random < 0.5)
            random = -1.0 + random;
	printf("THE random number is %f\n", random);
        Rbeta = random * factor;

        random = ( (double) rand() / RAND_MAX) ;
        if (random < 0.5)
            random = -1.0 + random;
	printf("THE random number is %f\n", random);
        Rgamma = random;

        random = ( (double) rand() / RAND_MAX) ;
        if (random < 0.5)
            random = -1.0 + random;
	printf("THE random number is %f\n", random);
        Rdelta = random * factor;

        for (i = 0; i < num; i++)
        for (j = 0; j < num; j++)
        {
    
            NEM2[i][j] = Ralpha * Ralpha * EM21[i][j] + Ralpha * Rgamma * EM1[i][j] + Rgamma * Rgamma * EM0[i][j];
        }
        for (i = 0; i < num; i++)
        {
            ipvt[i] = 0;
        }
        Icond = 0.0;
        det1[0] = det1[1] = 0.0;


        for (i = 0; i < num; i++)
        for (j = 0; j < num; j++)
	MAT_NUM[i*num+j] = NEM2[i][j];
        dgeco_(MAT_NUM,&dim,&dim,ipvt,&Icond,z);
        for (i = 0; i < num; i++)
        for (j = 0; j < num; j++)
	NEM2[i][j] = MAT_NUM[i*num+j];
        printf("Icond = %f \n",Icond);
        if (Icond >= EPS)
        {
        for (i = 0; i < num; i++)
        for (j = 0; j < num; j++)
	MAT_NUM[i*num+j] = NEM2[i][j];
            dgedi_(MAT_NUM,&dim,&dim,ipvt,det1,work,&job);
        for (i = 0; i < num; i++)
        for (j = 0; j < num; j++)
	NEM2[i][j] = MAT_NUM[i*num+j];
            for (i = 0; i < num; i++)
            for (j = 0; j < num; j++)
            {
                NEM1[i][j] = 2.0 * Ralpha * Rbeta * EM21[i][j] + (Ralpha * Rdelta + Rgamma * Rbeta) * EM1[i][j] + 2.0 * Rgamma * Rdelta * EM0[i][j];
                NEM0[i][j] = Rbeta * Rbeta * EM21[i][j] +  Rbeta * Rdelta * EM1[i][j] + Rdelta * Rdelta * EM0[i][j];
            }
           return(1);
        }
        count++;
        factor += 4.0;
    }
    return (0);
}
