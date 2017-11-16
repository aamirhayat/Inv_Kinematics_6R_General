/* 
This file contains the routine for reading the input parameter for the 
given robot. It consists of the fixes parameters like length of links,
offsets, twist angles and the RHS matrix for the inverse kinematics.


We use Fortran routine for reading, since there seems to be a problem in 
reading double precision input using a C compiler.
*/


#include <stdio.h>
#include "function.h"
#include <math.h>
#define SCALE 1.0

/*
    The global variables. Follow mostly to Raghavan/Roth paper.
*/

/*
     The L's are the cosines and M's are the sines of the twist angles.
*/
    extern coordinate L1,L2,L3,L4,L5,L6;
    extern coordinate M1,M2,M3,M4,M5,M6;

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
    extern coordinate Pi;

double rand1(seed) /* A modified Random Number Generator */
   int seed;
{
    srand(seed);
    return ((double) (rand() / 2147483647.0));
}

Read_Example()
{

    coordinate al1, al2, al3, al4, al5, al6;
    coordinate A1[4][4], A2[4][4], A3[4][4], A4[4][4], A5[4][4], A6[4][4];
    coordinate t1,t2,t3,t4,t5,t6;

    coordinate C1, C2, C3, C4, C5, C6;
    coordinate S1, S2, S3, S4, S5, S6;
    double Temp;

    coordinate sum, LHS[4][4], temp[4][4], temp1[4][4];

    int i,j,k,seed;
	FILE *fptr;

    if ((fptr = fopen("input.dat","r")) == NULL){
          printf("Error! opening file or No Input file found in ");
		// Program exits if the file pointer returns NULL.
          exit(1);
	} 
    fscanf(fptr,"%lf",&Temp);
    AA1 = (double) Temp * SCALE;
    fscanf(fptr,"%lf",&Temp);
    AA2 = (double) Temp * SCALE;
    fscanf(fptr,"%lf",&Temp);
    AA3 = (double) Temp * SCALE;
    fscanf(fptr,"%lf",&Temp);
    AA4 = (double) Temp * SCALE;
    fscanf(fptr,"%lf",&Temp);
    AA5 = (double) Temp * SCALE;
    fscanf(fptr,"%lf",&Temp);
    AA6 = (double) Temp * SCALE;
 
    fscanf(fptr,"%lf",&Temp);
    D_d1 = (double) Temp * SCALE;
    fscanf(fptr,"%lf",&Temp);
    D_d2 = (double) Temp * SCALE;
    fscanf(fptr,"%lf",&Temp);
    D_d3 = (double) Temp * SCALE;
    fscanf(fptr,"%lf",&Temp);
    D_d4 = (double) Temp * SCALE;
    fscanf(fptr,"%lf",&Temp);
    D_d5 = (double) Temp * SCALE;
    fscanf(fptr,"%lf",&Temp);
    D_d6 = (double) Temp * SCALE;
 
/*
    read_a_(&AA1, &AA2, &AA3, &AA4, &AA5, &AA6);

    read_d_(&D_d1, &D_d2, &D_d3, &D_d4, &D_d5, &D_d6);

*/
    //fprintf(stderr,"%g %g %g %g %g %g \n",AA1,AA2,AA3,AA4,AA5,AA6);
    //fprintf(stderr,"%g %g %g %g %g %g \n",D_d1,D_d2,D_d3,D_d4,D_d5,D_d6);

    fscanf(fptr,"%lf",&Temp);
    al1 = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    al2 = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    al3 = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    al4 = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    al5 = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    al6 = (double) Temp;

/*
    read_al_(&al1, &al2, &al3, &al4, &al5, &al6);
*/

    L1 = cos(al1 * Pi/180.0);
    L2 = cos(al2 * Pi/180.0);
    L3 = cos(al3 * Pi/180.0);
    L4 = cos(al4 * Pi/180.0);
    L5 = cos(al5 * Pi/180.0);
    L6 = cos(al6 * Pi/180.0);

    M1 = sin(al1 * Pi/180.0);
    M2 = sin(al2 * Pi/180.0);
    M3 = sin(al3 * Pi/180.0);
    M4 = sin(al4 * Pi/180.0);
    M5 = sin(al5 * Pi/180.0);
    M6 = sin(al6 * Pi/180.0);


/*
  Lets just try a good many random values of the manipulator parameters.

    read_rhs_(&lx,&mx,&nx,&ly,&my,&ny,&lz,&mz,&nz,&qx,&qy,&qz);
*/
    fscanf(fptr,"%lf",&Temp);
    lx = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    mx = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    nx = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    ly = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    my = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    ny = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    lz = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    mz = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    nz = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    qx = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    qy = (double) Temp;
    fscanf(fptr,"%lf",&Temp);
    qz = (double) Temp;

/*
     Now we compute the RHS matrix.
*/


/*
   Check whether the RHS orientation matrix is orthogonal.
*/
   
   
/*
   printf("0,1,  %g \n", lx * mx + ly * my + lz * mz);
   printf("0,2,  %g \n", lx * nx + ly * ny + lz * nz);
   printf("1,1,  %g \n", mx * mx + my * my + mz * mz);
   printf("1,2,  %g \n", mx * nx + my * ny + mz * nz);
   printf("2,2,  %g \n", nx * nx + ny * ny + nz * nz);
*/

   fclose(fptr);
}

