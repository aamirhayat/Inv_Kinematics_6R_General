/*
   This file consists of routines for reducing the problem to an eigenvalue 
   problem. It also involves transformations on the matrix such that the 
   resulting  matrix obtained for eigenvalue transformation has a low 
   condition number.
*/

#include <stdio.h>
#include "function.h"
#include <math.h>
#include <stdlib.h> 

#define MAX_TRANS 3

/* used by procedure gauss */
#define MAXG 24 
#define MAXG2 48

#define MAXD  16
#define MAXD2  32
/* we treat a number as real if its imaginary is smaller than this*/
#define IMAGINARY 1.0e-3 

/*
    The global variables. Follow mostly to Raghavan/Roth paper.
*/

/*
     The L's are the cosines and M's are the sines of the twist angles.
*/
    extern coordinate L1,L2,L3,L4,L5,L6;
    extern coordinate M1,M2,M3,M4,M5,M6;
    extern coordinate S1,S2,S3,S4,S5,S6;
    extern coordinate C1,C2,C3,C4,C5,C6;
    extern coordinate Verify_Results();

/*
      a_i is the length of link i.

      d_i is the offset distance of joint i;
*/
     extern coordinate AA1,AA2,AA3,AA4,AA5,AA6;
     extern coordinate D_d1,D_d2,D_d3,D_d4,D_d5,D_d6;

/*
    These are the entries of the RHS matrix.

      [ lx  mx nx qx ]
      [ ly  my ny qy ]
      [ lz  mz nz qz ]
      [ 0   0  0  1  ]

*/
    extern coordinate lx,ly,lz,mx,my,mz,nx,ny,nz,qx,qy,qz;
    extern coordinate u,v,w,p,q,r;


    extern coordinate Pi; 
    extern coordinate Rcond;

    extern struct M_ENTRY IRL[2][9];
    extern struct M_ENTRY IRRL[6][9]; /* A conversion matrix used in back substituting */

/*
    extern struct M_ENTRY Lpivot[8][9];
    extern coordinate Rpivot[8][6];
    extern int Ipivot[6];
*/
    
/*
    These matrices are generated from FMAT. In the process of reducing the 
    problem to an Eigenvalue problem. 

    EM2 : Coefficient of the quadratic term.
    EM1 : Coefficient of the linear term.
    EM0 : Coefficient of the constant terms.

    IEM1 = EM2^-1 * EM1;
    IEM0 = EM2^-1 * EM0;
*/
    extern coordinate EM2[12][12],EM21[12][12];
    extern coordinate EM1[12][12], IEM1[12][12];
    extern coordinate EM0[12][12], IEM0[12][12];

/*
 These matrices are used in this file for performing transformations.
*/
    coordinate NEM2[12][12],NEM21[12][12];
    coordinate NEM1[12][12], NIEM1[12][12];
    coordinate NEM0[12][12], NIEM0[12][12];
    coordinate work[500];

/*
 Used for Condition Number Evaluations in LAPACK routines .
*/
   coordinate vl[24][24],VR[24][24],rconde[24],rcondv[24],scale[28],abnrm;
/*
   EIG is the final matrix whose eigenvalues we compute to solve for x3.
*/
    coordinate EIG[24][24], EIG1[24][24],EIGcopy[24][24],EIG1copy[24][24];
    coordinate TT[24][24], TT1[24][24];
    coordinate SVD1[12][12],SVD[12][12];
/*
    int  fGeneralized;
   fGeneralized is a flag, set to 1, if Generalized Eigenvalue problem is
   used.
*/
  /* They contain the eigenvalues and   Eigenvectors. */

    extern coordinate D_wr[24],D_wi[24],D_beta[24],D_zz[24][24];
    extern coordinate Solution[6][16];
    extern int num_sols;
    extern int fGeneralized;

/* 
  These are the random elements chosen in the matrix transformations.
*/

    coordinate Ralpha, Rbeta, Rgamma, Rdelta;
    coordinate Icond; /* the condition number in the computations. */

/* 
  Mat_Mul: to multiply two matrices.
*/


Mat_MulG(A,B,C,order)
coordinate A[MAXG][MAXG];
coordinate B[MAXG][MAXG];
coordinate C[MAXG][MAXG];
int order;
{
    int i,j,k;
    coordinate sum;

    for (i = 0; i < order; i++)
        for (j = 0; j < order; j++)
        {
            sum = 0.0;
            for (k = 0; k < order; k++)
                sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
}


Mat_Mul(A,B,C,order)
coordinate A[12][12];
coordinate B[12][12];
coordinate C[12][12];
int order;
{
    int i,j,k;
    coordinate sum;

    for (i = 0; i < order; i++)
        for (j = 0; j < order; j++)
        {
            sum = 0.0;
            for (k = 0; k < order; k++)
                sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
}


Mat_MulD(A,B,C,order)
coordinate A[MAXD][MAXD];
coordinate B[MAXD][MAXD];
coordinate C[MAXD][MAXD];
int order;
{
    int i,j,k;
    coordinate sum;

    for (i = 0; i < order; i++)
        for (j = 0; j < order; j++)
        {
            sum = 0.0;
            for (k = 0; k < order; k++)
                sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
}
/* Given the sines and cosines of the angles, it performs an accurate computation
   of the solutions.
*/

Compute_Solutions(sol_index)
int sol_index;
{
    int i;
    coordinate temp;

    temp = asin(S3);
    if (C3 < 0.0)
        if (temp > 0.0)
            temp = Pi - temp;
        else 
             temp = -Pi - temp;
     Solution[2][sol_index] = temp;

     temp = asin(S4);
     if (C4 < 0.0)
         if (temp > 0.0)
             temp = Pi - temp;
         else 
             temp = -Pi - temp;
      Solution[3][sol_index] = temp;

     temp = asin(S5);
     if (C5 < 0.0)
         if (temp > 0.0)
             temp = Pi - temp;
         else 
             temp = -Pi - temp;
      Solution[4][sol_index] = temp;

     temp = asin(S1);
     if (C1 < 0.0)
         if (temp > 0.0)
             temp = Pi - temp;
         else 
             temp = -Pi - temp;
      Solution[0][sol_index] = temp;

     temp = asin(S2);
     if (C2 < 0.0)
         if (temp > 0.0)
             temp = Pi - temp;
         else 
             temp = -Pi - temp;
      Solution[1][sol_index] = temp;

     temp = asin(S6);
     if (C6 < 0.0)
         if (temp > 0.0)
             temp = Pi - temp;
         else 
             temp = -Pi - temp;
      Solution[5][sol_index] = temp;
}


/*
  Perform_Trans(): Performs a transformation on the x3 such that the resulting
  matrix obtained has a low condition number for reduction to an Eigenvalue
  problem.
*/

int Perform_Trans()

{
    coordinate EPS = 0.1;
    int i,j,k;
    coordinate temp1, temp2;
    int seed;
    coordinate random, z[12];
    extern float urand_();
    int count, dim, ipvt[12], job;
    coordinate det1[2];
    coordinate work[12];
    coordinate factor;

    seed = 10;
    count = 0;
    dim = 12;
    job = 11;
    factor = 1.0;
    srand(seed);
    while ( count < MAX_TRANS)
    { 
        random = ( (double) rand() / RAND_MAX) ;
        if (random < 0.5)
            random = -1.0 + random;
        Ralpha = random;
/*
	printf("THE random number is %f\n", random);
*/

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

        for (i = 0; i < 12; i++)
        for (j = 0; j < 12; j++)
        {
    
            NEM2[i][j] = Ralpha * Ralpha * EM21[i][j] + Ralpha * Rgamma * EM1[i][j] + Rgamma * Rgamma * EM0[i][j];
        }
        for (i = 0; i < 12; i++)
        {
            ipvt[i] = 0;
        }
        Icond = 0.0;
        det1[0] = det1[1] = 0.0;
        dgeco_(NEM2,&dim,&dim,ipvt,&Icond,z);
        printf("Icond = %f \n",Icond);
        if (Icond >= EPS)
        {
            dgedi_(NEM2,&dim,&dim,ipvt,det1,work,&job);
            for (i = 0; i < 12; i++)
            for (j = 0; j < 12; j++)
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

power()
{
        int i,j,k,l,m,n;
        int dim,job;
        int ipvt[24];
        coordinate det1[2];
        coordinate PA[24][24],PB[24][24],PX[24][24],test[24][24];
        coordinate vector[24],vectorcopy[24],z[24];
        coordinate work[24];
        coordinate EPS,lamdba,perturb,Rcond;
        coordinate temp,scale,la,lb;
	int 	index,ITERATIONS;


	/* initialize variables */
        dim = 24;
        EPS = 0.2;
        det1[0]=0;det1[1]=0;
	ITERATIONS = 22;
	perturb =1;

        for (i=0;i<24;i++)
        ipvt[i]=0;

	if (!fGeneralized)
	   {
	       for (i=0;i<24;i++)
		  {
		      if (fabs(D_wi[i]/D_beta[i]) <= EPS)
			 {
			     lamdba = perturb*(D_wr[i]+D_wi[i]);
			     printf(" This is the perturbed value %f\n",lamdba);
			     for (j=0;j<24;j++)
				for (k=0;k<24;k++)
				   PA[k][j] = EIGcopy[j][k];
			     for (j=0;j<24;j++)
				PA[j][j]= EIGcopy[j][j] - lamdba;
				
			     dgeco_(PA,&dim,&dim,ipvt,&Rcond,z);
			     printf("rcond for inversion in power= %15.20f \n",Rcond);

			     for (k=0; k<ITERATIONS; k++)
				{
				    /* normalize the eigenvector */
				    temp=0;index=0;
				    /* printf(" Here is the eigenvectors \n");
				    for (l=0;l<24;l++)
				    printf(" %15.20f \n",D_zz[i][l]);
				    */
				    for (l=0;l<24;l++)
				       if (((D_zz[i][l]>temp) && ( temp >= 0)) 
					 || ((D_zz[i][l] < temp ) && (temp <= 0)))
							{
							temp= D_zz[i][l];
							index = l;
							}
				    for (l=0;l<24;l++)
				       D_zz[i][l] /= temp;
				   
				    /*printf(" Here is the eigenvectors \n");
				    for (l=0;l<24;l++)
				    printf(" %15.20f \n",D_zz[i][l]);
				    */
				    job=1;
				    dgesl_(PA,&dim,&dim,ipvt,D_zz[i],&job);
				} 
				    printf(" This is the old eigenvalue %.20f\n", D_wr[i]);
				    D_wr[i] += 1/D_zz[i][index]; 
				    printf(" This is the new eigenvalue %.20f\n", D_wr[i]);
			 }
		      
		  }
	 }

	else
	   {
	       for (i=0;i<24;i++)
		  {
		      if (fabs(D_wi[i]/D_beta[i]) <= EPS)
			 {
			    perturb = 0.9;
			     lamdba = perturb*D_wr[i]/D_beta[i];
			     for (j=0;j<24;j++)
				for (k=0;k<24;k++)
				   {
				       PA[k][j] = EIGcopy[j][k];
				       PB[k][j] = EIG1copy[j][k];
				     PX[k][j] = PA[k][j]- lamdba*EIG1copy[j][k];
				   }

			     dgeco_(PX,&dim,&dim,ipvt,&Rcond,z);
			     printf("rcond for inversion in power= %15.10f \n",Rcond);
			     
			     for (k=0;k<ITERATIONS;k++)
				{

				    for (l=0;l<24;l++)
				       if (D_zz[i][l]>temp) temp= D_zz[i][l];
				    for (l=0;l<24;l++)
				       D_zz[i][l] /= temp;

				    for (l=0;l<24;l++)
				       { vector[l] = 0;
					 for (j=0;j<24;j++)
					    vector[l] += PB[l][j]*D_zz[i][j];
				     }

				    job=1;
				    dgesl_(PX,&dim,&dim,ipvt,vector,&job);

				    for (l=0;l<24;l++)
				       D_zz[i][l] = vector[l];
				}
			     /* get the approximation of eigenvalue */
			     for ( k=0; k<24;k++)
			       {
				vector[k]= 0;
				for (l =0; l< 24; l++)
				vector[k] += D_zz[i][l]*PA[l][k]; 
				}
				la= 0;
			     for (k=0;k<24;k++)
			      la+= vector[k] * D_zz[i][k];
			     for ( k=0; k<24;k++)
			       {
				vector[k]= 0;
				for (l =0; l< 24; l++)
				vector[k] += D_zz[i][l]*PB[l][k]; 
				}
				lb= 0;
			     for (k=0;k<24;k++)
			      lb+= vector[k] * D_zz[i][k];
			      D_wr[i]= la/lb;
				printf(" This is the new eigenvalue %f\n", D_wr[i]);

			 }
		  }
	   }
	
}







/* 
  Given EM2, EM1, EM0, it reduces the univariate equation solving to an 
  eigenvalue problem.
*/

/* set_up_degen moved to set_up.c */

/* 
  Given EM2, EM1, EM0, it reduces the univariate equation solving to an 
  eigenvalue problem.
*/

/* reduce_Eigen() moved to reduce_eigen.c */

/*
  Check_Eigenvector analyzes the eigenvector for some numerical properties. 
In particular it looks for the coefficients of maximum magnitude, so as to
minimize the effect of numerical errors.
*/

Check_Eigenvector(eig,row,px4,px5)
coordinate eig; /* the eigenvalue */
int row; /* its index, the eigenvector corresponding to this eigenvalue is D_zz[row] */
coordinate *px4, *px5;
{
    int i,j, col;
    coordinate  temp, temp1, max;

    col = 0;
  /*  if (fabs(eig) > 1.0)
        col = 12; */

   /*   Find the max. */
 /*   temp = fabs(D_zz[row][col]);
    j = 0;
    for (i = 1; i < 12; i++) 
    {
        temp1 = fabs(D_zz[row][col+i]);
        if (temp1 > temp)
        {
            temp = temp1;
            j = i;
        }
    } */
    j=0;
      
    switch (j)
    {

        case 0 : 
        
            *px5 = D_zz[row][col]/D_zz[row][col+1];
            *px4 = D_zz[row][col]/D_zz[row][col+3];
            return;

        case 1 : 
        
            *px5 = D_zz[row][col+1]/D_zz[row][col+2];
            *px4 = D_zz[row][col+1]/D_zz[row][col+4];
            return;

        case 2 : 
        
            *px5 = D_zz[row][col+1]/D_zz[row][col+2];
            *px4 = D_zz[row][col+2]/D_zz[row][col+5];
            return;

        case 3 : 
        
            *px5 = D_zz[row][col+3]/D_zz[row][col+4];
            *px4 = D_zz[row][col+3]/D_zz[row][col+6];
            return;

        case 4 : 
        
            *px5 = D_zz[row][col+4]/D_zz[row][col+5];
            *px4 = D_zz[row][col+4]/D_zz[row][col+7];
            return;

        case 5 : 
        
            *px5 = D_zz[row][col+4]/D_zz[row][col+5];
            *px4 = D_zz[row][col+5]/D_zz[row][col+8];
            return;

        case 6 : 
        
            *px5 = D_zz[row][col+6]/D_zz[row][col+7];
            *px4 = D_zz[row][col+6]/D_zz[row][col+9];
            return;

        case 7 : 
        
            *px5 = D_zz[row][col+7]/D_zz[row][col+8];
            *px4 = D_zz[row][col+7]/D_zz[row][col+10];
            return;

        case 8 : 
        
            *px5 = D_zz[row][col+8]/D_zz[row][col+11];
            *px4 = D_zz[row][col+10]/D_zz[row][col+11];
            return;

        case 9 : 
        
            *px5 = D_zz[row][col+9]/D_zz[row][col+10];
            *px4 = D_zz[row][col+6]/D_zz[row][col+9];
            return;

        case 10 : 
        
            *px5 = D_zz[row][col+10]/D_zz[row][col+11];
            *px4 = D_zz[row][col+7]/D_zz[row][col+10];
            return;

        case 11 : 
        
            *px5 = D_zz[row][col+10]/D_zz[row][col+11];
            *px4 = D_zz[row][col+8]/D_zz[row][col+11];
            return;

       default:   printf("j = %d \n",j);
   
       }

     /*  printf(" This is the eigenvectors \n");
       printf("  X4 = %f, X5 = %f\n",px4,px5); */
}

Compute_all_Solutions()

{
    int i,j,k,i1;
    coordinate EPS;
    coordinate temp,temp1,temp2,temp3;
    coordinate vector[9];
    coordinate lhs[6], rhs[6],tempv[6];
    coordinate rhs1, rhs2;
    coordinate c11,c12,c21,c22; 
    /* To solve the linear system of 2*2 equations. */ 
    coordinate eps_trig = 0.00005;
    coordinate det,x;

    if (Rcond < 0.0001)
        EPS = 0.000001 / Rcond;
    else
        EPS = 0.00005;
    EPS = 0.20005;
    num_sols = 0;
    /* power(); */
    for (i = 0; i < 24; i++)
        if (fabs(D_wi[i]/D_beta[i]) <= EPS)
        {
            if (fGeneralized)
            {
                temp1 = D_wr[i];
                temp2 = D_beta[i];
                S3 = 2 * temp1 * temp2 / ( temp2 * temp2 + temp1 * temp1);
                C3 = (temp2*temp2 - temp1 * temp1) / ( temp2 * temp2 + temp1 * temp1);
                Check_Eigenvector(temp1/D_beta[i],i,&temp2,&temp3);
            }
            else
            {
                temp1 = D_wr[i];
                S3 = 2 * temp1 / ( 1 + temp1 * temp1);
                C3 = (1 - temp1 * temp1) / ( 1 + temp1 * temp1);
                Check_Eigenvector(temp1,i,&temp2,&temp3);
            }

            S4 = 2 * temp2 / (1 + temp2 * temp2);
            C4 = (1 - temp2 * temp2) / (1 + temp2 * temp2);
            if ( fabs(1.0 - S4 * S4 - C4 * C4) > eps_trig)
                if (C4 == -1.0)
                    S4 = 0.0;
                else 
                {
                    x = S4 / (1 + C4);
                    S4 = 2 * x / (1 + x*x);
                    C4 = (1 - x*x) / (1 + x*x);
                }
            S5 = 2 * temp3 / (1 + temp3 * temp3);
            C5 = (1 - temp3 * temp3) / (1 + temp3 * temp3);
            if ( fabs(1.0 - S5 * S5 - C5 * C5) > eps_trig)
                if (C5 == -1.0)
                    S5 = 0.0;
                else 
                {
                    x = S5 / (1 + C5);
                    S5 = 2 * x / (1 + x*x);
                    C5 = (1 - x*x) / (1 + x*x);
                }

/*
   Compute the RHS vector entry for computing S1,S2, C1, C2.
*/
            vector[0] = S4 * S5;
            vector[1] = S4 * C5;
            vector[2] = C4 * S5;
            vector[3] = C4 * C5;
            vector[4] = S4;
            vector[5] = C4;
            vector[6] = S5;
            vector[7] = C5;
            vector[8] = 1;
            S1 = 0.0;
            C1 = 0.0;
            S2 = 0.0;
            C2 = 0.0;
            for (j = 0; j < 9; j++)
            {
                S1 = S1 + (IRL[0][j].Cos * C3 + IRL[0][j].Sin * S3 + IRL[0][j].Const) * vector[j];
                C1 = C1 + (IRL[1][j].Cos * C3 + IRL[1][j].Sin * S3 + IRL[1][j].Const) * vector[j];
            }

	     /*  printf(" This is vector number S1, %f\n",S1);
	       printf(" This is vector number C1 , %f\n",C1); */

            if ( fabs(1.0 - S1 * S1 - C1 * C1) > eps_trig)
                if (C1 == -1.0)
                    S1 = 0.0;
                else 
                {
                    x = S1 / (1 + C1);
                    S1 = 2 * x / (1 + x*x);
                    C1 = (1 - x*x) / (1 + x*x);
                }
	       /* printf(" This is improved vector number S1, %f\n",S1);
	       printf(" This is improved vector number C1 , %f\n",C1); */
/*

 Solve the upper triangular system for computing the RHS
*/
            for (j = 0; j < 9; j++)
            {
                S2 = S2 + (IRRL[4][j].Cos * C3 + IRRL[4][j].Sin * S3 + IRRL[4][j].Const) * vector[j];
                C2 = C2 + (IRRL[5][j].Cos * C3 + IRRL[5][j].Sin * S3 + IRRL[5][j].Const) * vector[j];
            }

	      /* printf(" This is vector number S2, %f\n",S2);
	       printf(" This is vector number C2 , %f\n",C2); */

	    for (i1=0;i1<6;i1++)
	     { 
	      tempv[i1] =0;
              for (j = 0; j < 9; j++)
	       tempv[i1] = tempv[i1] + (IRRL[i1][j].Cos * C3 + IRRL[i1][j].Sin * S3 + IRRL[i1][j].Const) * vector[j];

	      /* printf(" This is vector number [%d], %f\n",i1, tempv[i1]); */
                }

            if ( fabs(1.0 - S2 * S2 - C2 * C2) > eps_trig)
                if (C2 == -1.0)
                    S2 = 0.0;
                else 
                {
                    x = S2 / (1 + C2);
                    S2 = 2 * x / (1 + x*x);
                    C2 = (1 - x*x) / (1 + x*x);
                }
	     /*  printf(" This is improved vector number S2, %f\n",S2);
	       printf(" This is imporved vector number C2 , %f\n",C2);
	*/

/*
  The following code computes x6 from a 2 * 2 linear system. 
*/

            rhs1 = M3 * (S4 * C5 + S5 * C4 * L4) + S5 * L3 * M4;
            rhs2 = L5 * M3 *  (-S4 * S5 + C5 * C4 * L4) + L3 * (C5 * L5 * M4 + M5 * L4) - M3 * M4 * M5 * C4;
            c11 =  C2 *M2*M1*mz*L6 - S2*M2*C1*mx*L6 + S2*M2*C1*nx*M6 - L2*S1*M1*mx*L6- S2*M2*S1*my*L6+S2*M2*S1*ny*M6+L2*C1*M1*my*L6 - C2*M2*S1*L1*mx*L6 + C2*M2*S1*L1*nx*M6 - L2*L1*mz*L6 + C2*M2*C1*L1*my*L6 - C2*M2*C1*L1*ny*M6 - L2*C1*M1*ny*M6 + L2*S1*M1*nx*M6 - C2*M2*M1*nz*M6 + L2*L1*nz*M6;

            c12 = -C2*M2*M1*lz + S2*M2*C1*lx - L2*C1*M1*ly + L2*S1*M1*lx + S2*M2*S1*ly + L2*L1*lz - C2*M2*C1*L1*ly + C2*M2*S1*L1*lx;

            c21 =  -C2*M2*M1*lz + S2*M2*C1*lx - L2*C1*M1*ly + L2*S1*M1*lx + S2*M2*S1*ly + L2*L1*lz - C2*M2*C1*L1*ly + C2*M2*S1*L1*lx;
           
            c22 =  C2*M2*S1*L1*mx*L6 + S2*M2*C1*mx*L6 - C2*M2*C1*L1*my*L6 + C2*M2*C1*L1*ny*M6 + S2*M2*S1*my*L6 - C2*M2*M1*mz*L6 + C2*M2*M1*nz*M6 - S2*M2*C1*nx*M6 + L2*S1*M1*mx*L6 - L2*S1*M1*nx*M6 - S2*M2*S1*ny*M6 - L2*C1*M1*my*L6 + L2*C1*M1*ny*M6 - C2*M2*S1*L1*nx*M6 + L2*L1*mz*L6 - L2*L1*nz*M6;

            det = c11 * c22 - c12 * c21;
            if (det == 0.0)
                printf("Singular Matrix encountered in the computation of x6 \n");
            C6 = (c11 * rhs2 - c21 * rhs1) / det;
            S6 = (c22 * rhs1 - rhs2 * c12) / det;
            if ( fabs(1.0 - S6 * S6 - C6 * C6) > eps_trig)
                if (C6 == -1.0)
                    S6 = 0.0;
                else 
                {
                    x = S6 / (1 + C6);
                    S6 = 2 * x / (1 + x*x);
                    C6 = (1 - x*x) / (1 + x*x);
                }
            Verify_Results();
                   Improve_solutions(); 
/* Using Newton's method */
            Compute_Solutions(num_sols);
            num_sols++;
	    }
}


Improve_Eigen(Mat,B,lam,Vec)
coordinate Mat[24][24],B[24][24],lam,Vec[24];

{

    int i,j,k;
    int count, dim, ipvt[24],job,lda;
    coordinate A[24][24], rcond,z[24],v[24],w[24],sum;

    lda = 24;

    for (i = 0; i < 24; i++)
       for (j = 0; j < 24; j++)
           A[i][j] = Mat[i][j];

    dgeco_(A,&lda,&lda,ipvt,&rcond,z);
    printf(" IMPROVE: rcond = %12.8f \n",rcond);

    for (k = 0; k < 10; k++)
    {
        for (i = 0; i < 24; i++)
        {
            sum = 0.0;
            for (j = 0; j < 24; j++)
                sum = sum + B[i][j] * Vec[j];
            v[i] = sum;
        }

        dgesl_(A,&lda,&lda,ipvt,v,&job);
        sum = 0.0;
        
        for (i = 0; i < 24; i++)
            sum = sum + v[i] * v[i];

        sum = sqrt(sum);

        for (i = 0; i < 24; i++)
           Vec[i] = v[i]/sum;

        for(i = 0; i < 24; i++)
        {
            sum = 0.0;
            for (j = 0; j < 24; j++)
                sum = sum + Mat[i][j] * Vec[j];
            w[i] = sum;
        }

        lam = 0.0;
        for (i = 0; i < 24; i++)
            lam = lam + w[i] * Vec[i];
       
        printf("LAM = %12.8f \n",lam);
    }
     
}
