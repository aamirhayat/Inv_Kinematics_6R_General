/*
   This file consists of routines for reducing the problem to an eigenvalue 
   problem. It also involves transformations on the matrix such that the 
   resulting  matrix obtained for eigenvalue transformation has a low 
   condition number.
*/


/* used by procedure gauss */
/* we treat a number as real if its imaginary is smaller than this*/

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
