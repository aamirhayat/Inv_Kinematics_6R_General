
/*
    This file contains routines for improving the accuracy of the 
    results using Newton's method.
*/
#include <stdio.h>
#include "function.h"

    
    extern coordinate L1,L2,L3,L4,L5,L6;
    extern coordinate M1,M2,M3,M4,M5,M6;
    extern coordinate S1,S2,S3,S4,S5,S6;
    extern coordinate C1,C2,C3,C4,C5,C6;

    extern coordinate AA1,AA2,AA3,AA4,AA5,AA6;
    extern coordinate D_d1,D_d2,D_d3,D_d4,D_d5,D_d6;

    extern coordinate lx,ly,lz,mx,my,mz,nx,ny,nz,qx,qy,qz;
    extern coordinate u,v,w,p,q,r;

/*
  It is assumed that we are all the time improving on  values defined as tan x/2.
*/

Evaluate_Function(X,Vec,n)
coordinate *X, *Vec;  /* the vector X contains the set of values at which the 
                         function needs to be evaluated and Vec has the answer.
                       */
int n; /* The number of elements */

{
    coordinate n1,n2,n3,h1,h2,h3,m1,m2,m3,g1,g2,g3,r1,r2,r3,f1,f2,f3;
    coordinate S1,S2,S3,S4,S5, C1,C2,C3,C4,C5;

    coordinate x1,x2,x3,x4,x5;
    coordinate temp, sum;

    x1 = X[0];
    x2 = X[1]; 
    x3 = X[2];
    x4 = X[3];
    x5 = X[4];

/*
    printf(" IN COMPUTE_FUNCTION:X \n %18.14g, %18.14g, %18.14g, %18.14g, %18.14g \n", x1,x2,x3,x4,x5);
*/
    sum = 1 + x1 * x1;
    S1 = 2 * x1 / sum;
    C1 = (1 - x1 * x1) / sum;

    sum = 1 + x2 * x2;
    S2 = 2 * x2 / sum;
    C2 = (1 - x2 * x2) / sum;

    sum = 1 + x3 * x3;
    S3 = 2 * x3 / sum;
    C3 = (1 - x3 * x3) / sum;

    sum = 1 + x4 * x4;
    S4 = 2 * x4 / sum;
    C4 = (1 - x4 * x4) / sum;

    sum = 1 + x5 * x5;
    S5 = 2 * x5 / sum;
    C5 = (1 - x5 * x5) / sum;

    n1 = C1 * u + S1 * v;
    n2 = -L1 * (S1 * u - C1 * v) + M1 * w;
    n3 = M1 * (S1 * u - C1 * v) + L1 * w;

    h1 = C1 * p + S1 * q - AA1;
    h2 = - L1 * (S1 * p - C1 * q) + M1 * (r - D_d1);
    h3 = M1 * (S1 * p - C1 * q) + L1 * (r - D_d1);

    m1 = S5 * M5;
    m2 = C5 * L4 * M5 + M4 * L5;
    m3 = -C5 * M4 * M5 + L4 * L5;

    g1 = C5 * AA5 + AA4;
    g2 = -S5 * L4 * AA5 + M4 * D_d5;
    g3 = S5 * M4 * AA5 + L4 * D_d5 + D_d4;

    r1 = C4 * m1 + S4 * m2;
    r2 = -L3 * (S4 * m1 - C4 * m2) + M3 * m3;
    r3 = M3 * (S4 * m1 - C4 * m2) + L3 * m3;

    f1 = C4 * g1 + S4 * g2 + AA3;
    f2 = - L3 * (S4 * g1 - C4 * g2) + M3 * g3;
    f3 = M3 * (S4 * g1 - C4 * g2) + L3 * g3 + D_d3;

    Vec[0] = C3 * f1 + S3 * f2 - C2 * h1 - S2 * h2 + AA2;
    Vec[1] = S3 * f1 - C3 * f2 + L2 * (S2 * h1 - C2 * h2) - M2 * (h3 - D_d2);
    Vec[2] = f3 - M2 * (S2 * h1 - C2 * h2) - L2 * (h3 - D_d2);
    Vec[3] =  C3 * r1 + S3 * r2 - C2 * n1 - S2 * n2;
    Vec[4] =  S3 * r1 - C3 * r2 + L2 * (S2 * n1 - C2 * n2) - M2 * n3;
    x5 =  r3 - M2 * (S2 * n1 - C2 * n2) - L2 * n3;

/*
    printf(" IN COMPUTE_FUNCTION \n\n %18.14g, %18.14g, %18.14g, %18.14g, %18.14g, %18.14g \n", Vec[0], Vec[1], Vec[2], Vec[3], Vec[4],x5);
*/
}

/*  Given the set of values, it evaluates the Jacobian of the given system 
    of equations. */ 

Evaluate_Jacobian(X,Jac,n)
coordinate *X, Jac[5][5];/* the vector X contains the set of values at which 
                            the function needs to be evaluated and Vec has the
                            answer.*/
int n; /* The number of elements */


{
    int i,j,k;
/*
 These variable follow the same notation as in Raghavan/Roth paper on page 315. */

    coordinate n1,n2,n3,h1,h2,h3,m1,m2,m3,g1,g2,g3,r1,r2,r3,f1,f2,f3;
    coordinate S1,S2,S3,S4,S5, C1,C2,C3,C4,C5;

    coordinate x1,x2,x3,x4,x5;
    coordinate x12,x22,x32,x42,x52;
    coordinate sum1,sum2,sum3,sum4,sum5;
    coordinate temp, sum;

    coordinate h1x1, h2x1, h3x1; /* h1x1  is partial of h1 w.r.t. x1 and so on ...*/
    coordinate n1x1, n2x1, n3x1; /* n1x1  is partial of n1 w.r.t. x1 and so on ...*/

    coordinate m1x5, m2x5, m3x5; /* m1x5  is partial of m1 w.r.t. x5 and so on ...*/
    coordinate g1x5, g2x5, g3x5; /* g1x5  is partial of g1 w.r.t. x5 and so on ...*/

    coordinate r1x4,r1x5,r2x4,r2x5,r3x4,r3x5;
    coordinate f1x4,f1x5,f2x4,f2x5,f3x4,f3x5;

    x1 = X[0];
    x2 = X[1]; 
    x3 = X[2];
    x4 = X[3];
    x5 = X[4];

    x12 = x1 * x1;
    x22 = x2 * x2;
    x32 = x3 * x3;
    x42 = x4 * x4;
    x52 = x5 * x5;

    sum1 = 1 + x12;
    S1 = 2 * x1 / sum1;
    C1 = (1 - x12) / sum1;

    sum2 = 1 + x22;
    S2 = 2 * x2 / sum2;
    C2 = (1 - x22) / sum2;

    sum3 = 1 + x32;
    S3 = 2 * x3 / sum3;
    C3 = (1 - x32) / sum3;

    sum4 = 1 + x42;
    S4 = 2 * x4 / sum4;
    C4 = (1 - x42) / sum4;

    sum5 = 1 + x52;
    S5 = 2 * x5 / sum5;
    C5 = (1 - x52) / sum5;

    n1 = C1 * u + S1 * v;
    n2 = -L1 * (S1 * u - C1 * v) + M1 * w;
    n3 = M1 * (S1 * u - C1 * v) + L1 * w;

    h1 = C1 * p + S1 * q - AA1;
    h2 = - L1 * (S1 * p - C1 * q) + M1 * (r - D_d1 );
    h3 = M1 * (S1 * p - C1 * q) + L1 * (r - D_d1 );

    m1 = S5 * M5;
    m2 = C5 * L4 * M5 + M4 * L5;
    m3 = -C5 * M4 * M5 + L4 * L5;

    g1 = C5 * AA5 + AA4;
    g2 = -S5 * L4 * AA5 + M4 * D_d5;
    g3 = S5 * M4 * AA5 + L4 * D_d5 + D_d4;

    r1 = C4 * m1 + S4 * m2;
    r2 = -L3 * (S4 * m1 - C4 * m2) + M3 * m3;
    r3 = M3 * (S4 * m1 - C4 * m2) + L3 * m3;

    f1 = C4 * g1 + S4 * g2 + AA3;
    f2 = - L3 * (S4 * g1 - C4 * g2) + M3 * g3;
    f3 = M3 * (S4 * g1 - C4 * g2) + L3 * g3 + D_d3;


    h1x1 =(-4 * x1 * p + 2 * (1 - x12) * q) / (sum1 * sum1); 
    h2x1 = -L1 * (4 * x1 * q + 2 * (1 - x12) * p) / (sum1 * sum1); 
    h3x1 = M1 * (4 * x1 * q + 2 * (1 - x12) * p) / (sum1 * sum1); 

    n1x1 =(-4 * x1 * u + 2 * (1 - x12) * v) / (sum1 * sum1); 
    n2x1 = -L1 * (4 * x1 * v + 2 * (1 - x12) * u) / (sum1 * sum1); 
    n3x1 = M1 * (4 * x1 * v + 2 * (1 - x12) * u) / (sum1 * sum1); 


    m1x5 = M5 * 2 * (1 - x52) / (sum5 * sum5);
    m2x5 = L4 * M5 * -4 * x5 / (sum5 * sum5);
    m3x5 = M4 * M5 * 4 * x5 / (sum5 * sum5);


    g1x5 = AA5 * -4 * x5 / (sum5 * sum5);
    g2x5 = -L4 * AA5 * 2 * (1 - x52) / (sum5 * sum5); 
    g3x5 =  M4 * AA5 * 2 * (1 - x52) / (sum5 * sum5); 


    r1x4 = (m1 * -4 * x4 + 2 * m2 * (1 - x42)) / (sum4 * sum4);
    r1x5 = C4 * m1x5 + S4 * m2x5;

    r2x4 = -L3 * (m2 * 4 * x4 + 2 * m1 * (1 - x42)) / (sum4 * sum4);
    r2x5 = -L3 * (S4 * m1x5 - C4 * m2x5) + M3 * m3x5;;

    r3x4 = M3 * (m2 * 4 * x4 + 2 * m1 * (1 - x42)) / (sum4 * sum4);
    r3x5 = M3 * (S4 * m1x5 - C4 * m2x5) + L3 * m3x5;;

    f1x4 = (-4 * x4 * g1 + 2 * (1 - x42) * g2) / (sum4 * sum4);
    f1x5 = C4 * g1x5 + S4 * g2x5;


    f2x4 = -L3 * (g2 * 4 * x4 + 2 * g1 * (1 - x42)) / (sum4 * sum4);
    f2x5 = -L3 * (S4 * g1x5 - C4 * g2x5) + M3 * g3x5;;

    f3x4 = M3 * (g2 * 4 * x4 + 2 * g1 * (1 - x42)) / (sum4 * sum4);
    f3x5 = M3 * (S4 * g1x5 - C4 * g2x5) + L3 * g3x5;;

/*
  We have to realize that the main routine for Newton's method is implemented
  in Fortran. As a result we will be entering the elements of the Jacobian 
  in column order format. */
     
  Jac[0][0] = - C2 * h1x1 - S2 * h2x1;
  Jac[1][0] = (h1 * 4 * x2 - h2 * 2 * (1 - x22)) / (sum2 * sum2);
  Jac[2][0] = (f1 * -4 * x3 + f2 * 2 * (1 - x32)) / (sum3 * sum3);
  Jac[3][0] = C3 * f1x4 + S3 * f2x4;
  Jac[4][0] = C3 * f1x5 + S3 * f2x5;


  Jac[0][1] = L2 * (S2 * h1x1 - C2 * h2x1) - M2 * h3x1;
  Jac[1][1] = L2 * (h1 * 2 * (1 - x22) + 4 * x2 * h2) / (sum2 * sum2);
  Jac[2][1] = (f1 * 2 * (1 - x32) + 4 * x3  * f2) / (sum3 * sum3);
  Jac[3][1] = S3  * f1x4 - C3 * f2x4;
  Jac[4][1] = S3 * f1x5 - C3 * f2x5;

  Jac[0][2] = -M2 * (S2 * h1x1 - C2 * h2x1) - L2 * h3x1;
  Jac[1][2] = -M2 * (h1 * 2 * (1 - x22) + 4 * x2 * h2) / (sum2 * sum2);
  Jac[2][2] = 0.0;
  Jac[3][2] = f3x4;
  Jac[4][2] = f3x5;

  Jac[0][3] = -1* (C2 * n1x1 + S2 * n2x1);
  Jac[1][3] = (n1 * 4 * x2 - n2 * 2 * (1 - x22)) / (sum2 * sum2);
  Jac[2][3] = (-4 * x3 * r1 + 2 * r2 * (1 - x32)) / (sum3 * sum3);
  Jac[3][3] = C3 * r1x4 + S3 * r2x4;
  Jac[4][3] = C3 * r1x5 + S3 * r2x5;


  Jac[0][4] = L2 * (S2 * n1x1 - C2 * n2x1) - M2 * n3x1;
  Jac[1][4] = L2 * (n1 * 2 * (1 - x22) + 4 * x2 * n2) / (sum2 * sum2);
  Jac[2][4] = (r1 * 2 * (1 - x32) + 4 * x3  * r2) / (sum3 * sum3);
  Jac[3][4] = S3  * r1x4 - C3 * r2x4;
  Jac[4][4] = S3 * r1x5 - C3 * r2x5;
}

/*
  This function: fcn: is called by hybrj for evaluating the given function and
  its Jacobians.
*/

fcn(pn,X,fvec,fjac,pldfjac,piflag)
int *pn,*pldfjac,*piflag;
coordinate X[5], fvec[5],fjac[5][5];

{
    int n,ldfjac,iflag;

    n = *pn;
    ldfjac = *pldfjac;
    iflag = *piflag;

    if (iflag == 1)
    {
        Evaluate_Function(X,fvec,n);
        return;
    }
    if (iflag == 2)
    {
        Evaluate_Jacobian(X,fjac,n);
        return;
    }
}

    
    

/*
   Improve_solutions(): Main Procedure for Newton's iteration based improvement.
*/

Improve_solutions()
{

    int n,ldfjac,maxfev,mode,nprint,info,nfev,njev,lr;
    coordinate xtol, factor;
    coordinate X[5],fvec[5],fjac[5][5],diag[5],r[20],qtf[5],wal[5],wAA1[5],wAA2[5],wAA3[5],wAA4[5];
    int *func;
    coordinate x1,x2,x3,x4,x5;
    
    n = 5;
    ldfjac = 5;

/* Compute the initial estimates */

    if (C1 == -1.0)
        X[0] = 1.0e16;
    else
        X[0] = S1 / (1 + C1);

    if (C2 == -1.0)
        X[1] = 1.0e16;
    else
        X[1] = S2 / (1 + C2);

    if (C3 == -1.0)
        X[2] = 1.0e16;
    else
        X[2] = S3 / (1 + C3);

    if (C4 == -1.0)
        X[3] = 1.0e16;
    else
        X[3] = S4 / (1 + C4);

    if (C5 == -1.0)
        X[4] = 1.0e16;
    else
        X[4] = S5 / (1 + C5);


    xtol = 0.0001;
    maxfev = 10;
    mode = 1;
    factor = 80.0;
    nprint = 0;
    lr = 20;

    hybrj_(fcn,&n,X,fvec,fjac,&ldfjac,&xtol,&maxfev,diag,&mode,&factor,&nprint, &info,&nfev,&njev,r,&lr,qtf,wAA1,wAA2,wAA3,wAA4);

    x1 = X[0];   
    x2 = X[1];   
    x3 = X[2];   
    x4 = X[3];   
    x5 = X[4];   
/*
    printf(" OUTPUT:X \n\n %18.14g, %18.14g, %18.14g, %18.14g, %18.14g \n", x1, x2, x3, x4, x5);
    printf(" OUTPUT:F \n\n %18.14g, %18.14g, %18.14g, %18.14g, %18.14g \n", fvec[0], fvec[1], fvec[2], fvec[3], fvec[4]);
*/

    S1 = 2 * x1 / (1 + x1 * x1);
    C1 = (1 - x1 * x1) / (1 + x1 * x1);
    S2 = 2 * x2 / (1 + x2 * x2);
    C2 = (1 - x2 * x2) / (1 + x2 * x2);
    S3 = 2 * x3 / (1 + x3 * x3);
    C3 = (1 - x3 * x3) / (1 + x3 * x3);
    S4 = 2 * x4 / (1 + x4 * x4);
    C4 = (1 - x4 * x4) / (1 + x4 * x4);
    S5 = 2 * x5 / (1 + x5 * x5);
    C5 = (1 - x5 * x5) / (1 + x5 * x5);
    if (info != 1)
        printf(" From Newton's Method: info = %d \n",info);
}
