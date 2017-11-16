#include <stdio.h>
#include "function.h"

/*
    The global variables. Follow notations from Raghavan/Roth paper.
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

/* 

   The u,v,w,p,q,r  correspond to the expressions defined in Raghavan/Roth 
   paper
*/
 
    extern coordinate u,v,w,p,q,r;


    extern coordinate Rsim1[6][2];
    extern coordinate Rsim2[8][6];
    extern struct M_ENTRY Lsim1[6][9];
    extern struct M_ENTRY Lsim2[8][9];


initialize_variables()
{


/*
        ai2 = ai * ai;
        di2 = di * di;

        D_d1r = D_d1 - r;
*/

    coordinate  AA12,AA22,AA32,AA42,AA52,AA62;
    coordinate  D_d12,D_d22,D_d32,D_d42,D_d52,D_d62;
    coordinate D_d1r,p2,q2,r2;;

    u = mx * M6 + nx * L6;
    v = my * M6 + ny * L6;
    w = mz * M6 + nz * L6;


    AA12 = AA1 * AA1;
    AA22 = AA2 * AA2;
    AA32 = AA3 * AA3;
    AA42 = AA4 * AA4;
    AA52 = AA5 * AA5;
    AA62 = AA6 * AA6;

    D_d12 = D_d1 * D_d1;
    D_d22 = D_d2 * D_d2;
    D_d32 = D_d3 * D_d3;
    D_d42 = D_d4 * D_d4;
    D_d52 = D_d5 * D_d5;
    D_d62 = D_d6 * D_d6;


    p = -lx * AA6 - u*D_d6 + qx;
    q = -ly * AA6 - v*D_d6 + qy;
    r = -lz * AA6 - w*D_d6 + qz;

    p2 = p * p;
    q2 = q * q;
    r2 = r * r;

    D_d1r = D_d1 - r;

    Rsim1[0][0] = M1 * p;
    Rsim1[0][1] = - M1 * q;
     
    Rsim1[1][0] = M1 * u;
    Rsim1[1][1] = - M1 * v;
     
    Rsim1[2][0] = -2 * q * AA1;
    Rsim1[2][1] = -2 * p * AA1 ;
     
    Rsim1[3][0] = -AA1 * v;
    Rsim1[3][1] = -AA1 * u;
     
    Rsim1[4][0] = (-q*w - D_d1r*v)*M1 - AA1*L1*u;
    Rsim1[4][1] = (-D_d1r*u - p*w)*M1 + AA1*L1*v;
     
    Rsim1[5][0] = (-2*q*v*p + r2*u+2*D_d1r*w*p - 2*r*D_d1*u + (AA12 + D_d12 - p2 + q2)*u)*M1 + (2*AA1*v*r - 2*q*AA1*w - 2*AA1*v*D_d1)*L1;
    Rsim1[5][1] = (-2*D_d1*w*q + 2*p*u*q - (AA12 - 2*r*D_d1)*v - (p2+r2)*v -(D_d12 - q2)*v + 2*r*w*q)*M1 + (-2*p*w - 2*u*D_d1 + 2*u*r)*AA1*L1;
     
    Rsim2[0][0] = -L1*p;
    Rsim2[0][1] = q;
    Rsim2[0][2] = L1*q;
    Rsim2[0][3] = p;
    Rsim2[0][4] = -M1*D_d1r;
    Rsim2[0][5] = -AA1;

    Rsim2[1][0] = q;
    Rsim2[1][1] = L1*p;
    Rsim2[1][2] = p;
    Rsim2[1][3] = -L1*q;
    Rsim2[1][4] = -AA1;
    Rsim2[1][5] = D_d1r*M1;

    Rsim2[2][0] = -L1*u;
    Rsim2[2][1] = v;
    Rsim2[2][2] = L1*v;
    Rsim2[2][3] = u;
    Rsim2[2][4] = M1*w;
    Rsim2[2][5] = 0.0;

    Rsim2[3][0] = v;
    Rsim2[3][1] = L1*u;
    Rsim2[3][2] = u;
    Rsim2[3][3] = -L1*v;
    Rsim2[3][4] = 0.0;
    Rsim2[3][5] = -M1*w;

    Rsim2[4][0] = (D_d1r*v + q*w)*L1 - AA1*M1*u;
    Rsim2[4][1] = D_d1r*u + p*w;
    Rsim2[4][2] = (D_d1r*u + p*w)*L1 + AA1*M1*v;
    Rsim2[4][3] = -q*w - D_d1r*v;
    Rsim2[4][4] = (q*u - p*v)*M1 - AA1*L1*w;
    Rsim2[4][5] = 0.0;


    Rsim2[5][0] = D_d1r*u + p*w;
    Rsim2[5][1] = (-D_d1r*v - q*w)*L1 + AA1*M1*u;
    Rsim2[5][2] = -D_d1r*v - q*w;
    Rsim2[5][3] = (-D_d1r*u - p*w)*L1 - AA1*M1*v;
    Rsim2[5][4] = 0.0;
    Rsim2[5][5] = (p*v - q*u)*M1 + AA1*L1*w;

    Rsim2[6][0] = (-q*w-v*D_d1r)*2*AA1*M1+((-r2+2*r*D_d1-q2-D_d12+p2-AA12)*u+2*(q*v+r*w-D_d1*w)*p)*L1;

    Rsim2[6][1] = 2*(D_d1r*w*q-r*D_d1*v)-2*p*u*q+(D_d12-q2+r2+p2-AA12)*v;

    Rsim2[6][2] = (-p*w-u*D_d1+u*r)*2*AA1*M1+(2*D_d1*w*q-2*r*w*q-2*p*u*q+(r2-q2-2*r*D_d1+D_d12+p2+AA12)*v)*L1;


    Rsim2[6][3] = -(2*r*D_d1+AA12)*u-2*(r*w+q*v-D_d1*w)*p+(-p2+r2+D_d12+q2)*u;

    Rsim2[6][4] = ((p2+AA12)*w+2*p*u*D_d1r+2*q*v*D_d1-(r2+D_d12-2*r*D_d1-q2)*w-2*q*v*r)*M1+(q*u-v*p)*2*AA1*L1;


    Rsim2[6][5] = 2*(q*v+p*u+r*w-D_d1*w)*AA1;


    Rsim2[7][0] = 2*(D_d1*w*q-r*w*q-r*D_d1*v)-2*p*u*q-(-D_d12+q2-r2-p2+AA12)*v;
    Rsim2[7][1] = 2*(q*w+v*D_d1-v*r)*AA1*M1-2*(q*v+r*w-D_d1*w)*p*L1+(-2*r*D_d1+r2+AA12+D_d12-p2+q2)*u*L1;

Rsim2[7][2] = (-2*r*D_d1-AA12)*u-2*(r*w+q*v-D_d1*w)*p+(-p2+r*r+D_d1*D_d1+q*q)*u;

    Rsim2[7][3] = (p*w-u*r+u*D_d1)*2*AA1*M1+(-2*(D_d1*w-p*u-r*w)*q+(-AA12+2*r*D_d1-p2-r2-D_d12+q2)*v)*L1;


    Rsim2[7][4] = 2*(q*v+p*u+r*w-D_d1*w)*AA1;

    Rsim2[7][5] = (2*(q*v*r-p*u*D_d1r-q*v*D_d1)+(-q2-AA12+D_d12-p2+r2-2*r*D_d1) *w)*M1+(v*p-q*u)*2*AA1*L1;

   

    Lsim1[0][0].Cos = 0.0;
    Lsim1[0][0].Sin = -M2*L4*AA5;
    Lsim1[0][0].Const = 0.0;

    Lsim1[0][1].Cos = M2*L3*AA5;
    Lsim1[0][1].Sin = 0.0;
    Lsim1[0][1].Const = L2*M3*AA5;

    Lsim1[0][2].Cos = M2*L3*L4*AA5;
    Lsim1[0][2].Sin = 0.0;
    Lsim1[0][2].Const = L2*M3*L4*AA5;

    Lsim1[0][3].Cos = 0.0;
    Lsim1[0][3].Sin = M2*AA5;
    Lsim1[0][3].Const = 0;


    Lsim1[0][4].Cos = AA4*L3*M2;
    Lsim1[0][4].Sin = M2*M4*D_d5;
    Lsim1[0][4].Const = L2*M3*AA4;

    Lsim1[0][5].Cos = -M2*L3*M4*D_d5;
    Lsim1[0][5].Sin = M2*AA4;
    Lsim1[0][5].Const = -L2*M3*M4*D_d5;

    Lsim1[0][6].Cos = -M2*M3*M4*AA5;
    Lsim1[0][6].Sin = 0.0;
    Lsim1[0][6].Const = L2*L3*M4*AA5;

    Lsim1[0][7].Cos = 0.0;
    Lsim1[0][7].Sin = 0.0;
    Lsim1[0][7].Const = 0.0;

    Lsim1[0][8].Cos = -M2*M3*D_d4-M2*M3*L4*D_d5;
    Lsim1[0][8].Sin = M2*AA3;
    Lsim1[0][8].Const = L1*D_d1-L1*r+L2*L3*L4*D_d5+L2*D_d3+L2*L3*D_d4+D_d2;



    Lsim1[1][0].Cos = M2*L3*M5;
    Lsim1[1][0].Sin = 0.0;
    Lsim1[1][0].Const = L2*M3*M5;

    Lsim1[1][1].Cos = 0.0;
    Lsim1[1][1].Sin = M2*L4*M5;
    Lsim1[1][1].Const = 0.0;

    Lsim1[1][2].Cos = 0.0;
    Lsim1[1][2].Sin = M2*M5;
    Lsim1[1][2].Const = 0.0;


    Lsim1[1][3].Cos = -M2*L3*L4*M5;
    Lsim1[1][3].Sin = 0.0;
    Lsim1[1][3].Const = -L2*M3*L4*M5;

    Lsim1[1][4].Cos = 0.0;
    Lsim1[1][4].Sin = M2*M4*L5;
    Lsim1[1][4].Const = 0.0;

    Lsim1[1][5].Cos = -M2*L3*M4*L5;
    Lsim1[1][5].Sin = 0.0;
    Lsim1[1][5].Const = -L2*M3*M4*L5;

    Lsim1[1][6].Cos = 0.0;
    Lsim1[1][6].Sin = 0.0;
    Lsim1[1][6].Const = 0.0;

    Lsim1[1][7].Cos = M2*M3*M4*M5;
    Lsim1[1][7].Sin = 0.0;
    Lsim1[1][7].Const = -L2*L3*M4*M5;

    Lsim1[1][8].Cos = -M2*M3*L4*L5;
    Lsim1[1][8].Sin = 0.0;  
    Lsim1[1][8].Const = L2*L3*L4*L5-L1*w;

    Lsim1[2][0].Cos = -2*L4*AA5*AA2;
    Lsim1[2][0].Sin = -2*M2*L4*AA5*D_d2;
    Lsim1[2][0].Const = -2*L4*AA5*AA3;

    Lsim1[2][1].Cos = 2*AA5*D_d2*L3*M2;
    Lsim1[2][1].Sin = -2*L3*AA5*AA2;
    Lsim1[2][1].Const = 2*M3*AA5*D_d3+2*L2*M3*AA5*D_d2;

    Lsim1[2][2].Cos = 2*M2*L3*L4*AA5*D_d2;
    Lsim1[2][2].Sin = -2*L3*L4*AA5*AA2;
    Lsim1[2][2].Const = 2*M3*L4*AA5*D_d3+2*L2*M3*L4*AA5*D_d2;

    Lsim1[2][3].Cos = 2*AA5*AA2;
    Lsim1[2][3].Sin = 2*M2*AA5*D_d2;
    Lsim1[2][3].Const = 2*AA5*AA3;


    Lsim1[2][4].Cos = 2*M4*D_d5*AA2+2*AA4*D_d2*L3*M2;
    Lsim1[2][4].Sin =  2*M2*M4*D_d5*D_d2-2*L3*AA4*AA2;
    Lsim1[2][4].Const = 2*M3*AA4*D_d3+2*L2*M3*AA4*D_d2+2*M4*D_d5*AA3;

    Lsim1[2][5].Cos = 2*AA4*AA2-2*M2*L3*M4*D_d5*D_d2;
    Lsim1[2][5].Sin =  2*M2*AA4*D_d2+2*L3*M4*D_d5*AA2;
    Lsim1[2][5].Const = -2*L2*M3*M4*D_d5*D_d2+2*AA4*AA3-2*M3*M4*D_d5*D_d3;

    Lsim1[2][6].Cos =  -2*M2*M3*M4*AA5*D_d2;
    Lsim1[2][6].Sin =  2*M3*M4*AA5*AA2;
    Lsim1[2][6].Const = 2*L3*M4*AA5*D_d3+2*M4*AA5*D_d4+2*L2*L3*M4*AA5*D_d2;

    Lsim1[2][7].Cos =  0.0;
    Lsim1[2][7].Sin = 0.0;
    Lsim1[2][7].Const = 2*AA5*AA4;

    Lsim1[2][8].Cos =  2*AA3*AA2-2*M2*M3*L4*D_d5*D_d2-2*M2*M3*D_d4*D_d2;
    Lsim1[2][8].Sin = 2*M3*L4*D_d5*AA2+2*M3*D_d4*AA2+2*M2*AA3*D_d2;
    Lsim1[2][8].Const = D_d3*D_d3-p*p-q*q-r*r-D_d1*D_d1-AA1*AA1+AA3*AA3+D_d4*D_d4+AA5*AA5+AA2*AA2+2*r*D_d1+D_d2*D_d2+AA4*AA4+D_d5*D_d5+2*L3*D_d4*D_d3+2*L3*L4*D_d5*D_d3+2*L2*L3*D_d4*D_d2+2*L4* D_d5*D_d4+2*L2*D_d3*D_d2+2*L2*L3*L4*D_d5*D_d2;



    Lsim1[3][0].Cos = D_d2*M2*L3*M5;
    Lsim1[3][0].Sin = -AA2*L3*M5;
    Lsim1[3][0].Const = D_d2*L2*M3*M5+D_d3*M3*M5;

    Lsim1[3][1].Cos = AA2*L4*M5;
    Lsim1[3][1].Sin = D_d2*M2*L4*M5;
    Lsim1[3][1].Const = AA3*L4*M5;

    Lsim1[3][2].Cos = AA2*M5;
    Lsim1[3][2].Sin = D_d2*M2*M5;
    Lsim1[3][2].Const = AA3*M5;


    Lsim1[3][3].Cos = -D_d2*M2*L3*L4*M5;
    Lsim1[3][3].Sin = AA2*L3*L4*M5;
    Lsim1[3][3].Const = -D_d3*M3*L4*M5-D_d2*L2*M3*L4*M5;

    Lsim1[3][4].Cos = AA2*M4*L5;
    Lsim1[3][4].Sin = D_d2*M2*M4*L5;
    Lsim1[3][4].Const = AA3*M4*L5;

    Lsim1[3][5].Cos = -D_d2*M2*L3*M4*L5;
    Lsim1[3][5].Sin = AA2*L3*M4*L5;
    Lsim1[3][5].Const = -D_d3*M3*M4*L5-D_d2*L2*M3*M4*L5;

    Lsim1[3][6].Cos = 0.0;
    Lsim1[3][6].Sin = 0.0;
    Lsim1[3][6].Const = AA4*M5;

    Lsim1[3][7].Cos = D_d2*M2*M3*M4*M5;
    Lsim1[3][7].Sin = -AA2*M3*M4*M5;
    Lsim1[3][7].Const = -D_d4*M4*M5-D_d2*L2*L3*M4*M5-D_d3*L3*M4*M5;

    Lsim1[3][8].Cos = -D_d2*M2*M3*L4*L5;
    Lsim1[3][8].Sin = AA2*M3*L4*L5;
    Lsim1[3][8].Const = D_d3*L3*L4*L5+D_d4*L4*L5+D_d2*L2*L3*L4*L5+D_d5*L5-p*u-r*w-q*v+D_d1*w;



    Lsim1[4][0].Cos = AA3*M2*M3*M5-AA2*L2*L3*M5-M2*L3*AA5*L5;
    Lsim1[4][0].Sin = L4*D_d5*M2*M5+D_d4*M2*M5+M2*D_d3*L3*M5;
    Lsim1[4][0].Const = -L2*AA3*L3*M5-M3*AA5*L2*L5+AA2*M2*M3*M5;

    Lsim1[4][1].Cos = -M2*L3*D_d4*L4*M5-M2*D_d3*L4*M5-D_d5*M2*L3*M5;
    Lsim1[4][1].Sin = -AA2*L2*L4*M5-M2*AA5*L4*L5+M2*AA4*M4*M5;
    Lsim1[4][1].Const = -L2*D_d5*M3*M5-M3*D_d4*L2*L4*M5;

    Lsim1[4][2].Cos = -M2*L3*L4*D_d5*M5-M2*D_d3*M5-M2*L3*D_d4*M5;
    Lsim1[4][2].Sin = -AA2*L2*M5-AA5*M2*L5;
    Lsim1[4][2].Const = -M3*L4*D_d5*L2*M5-M3*D_d4*L2*M5;


    Lsim1[4][3].Cos = M2*L4*(AA5*L3*L5-AA3*M3*M5)+(L3*M5)*(AA2*L2*L4-AA4*M2*M4);
    Lsim1[4][3].Sin = -D_d4*M2*L4*M5-M2*D_d5*M5-M2*D_d3*L3*L4*M5;
    Lsim1[4][3].Const = -L2*AA4*M3*M4*M5-AA2*M2*M3*L4*M5+L2*AA5*M3*L4*L5+L2*AA3*L3*L4*M5;

    Lsim1[4][4].Cos = -M2*L3*D_d4*M4*L5-M2*D_d3*M4*L5;
    Lsim1[4][4].Sin = -AA2*L2*M4*L5+M4*AA5*M2*M5-M2*AA4*L4*L5;
    Lsim1[4][4].Const = -M3*D_d4*L2*M4*L5;

    Lsim1[4][5].Cos = -M2*L3*M4*AA5*M5+AA4*M2*L3*L4*L5+AA2*L2*L3*M4*L5-AA3*M2*M3*M4*L5;
    Lsim1[4][5].Sin = -M2*D_d3*L3*M4*L5-D_d4*M2*M4*L5;
    Lsim1[4][5].Const = -M3*M4*AA5*L2*M5-AA2*M2*M3*M4*L5+L2* AA3*L3*M4*L5+L2*AA4*M3*L4*L5;

    Lsim1[4][6].Cos = M4*D_d5*M2*M3*M5;
    Lsim1[4][6].Sin = 0.0;
    Lsim1[4][6].Const = -L2*M4*D_d5*L3*M5;

    Lsim1[4][7].Cos = -AA3*M2*L3*M4*M5-M2*M3*AA5*M4*L5-M2* M3*AA4*L4*M5-AA2*L2*M3*M4*M5;
    Lsim1[4][7].Sin = M2*D_d3*M3*M4*M5;
    Lsim1[4][7].Const = -AA2*M2*L3*M4*M5+L3*AA5*L2*M4*L5+L3* AA4*L2*L4*M5-L2*AA3*M3*M4*M5;

    Lsim1[4][8].Cos = AA2*L2*M3*L4*L5-M2*M3*AA4*M4*L5-L4* AA5*M2*M3*M5+AA3*M2*L3*L4*L5;
    Lsim1[4][8].Sin = -M2*D_d3*M3*L4*L5;
    Lsim1[4][8].Const = L1*p*v-q*L1*u+L2*L4*AA5*L3*M5+AA2*M2*L3* L4*L5+L2*AA3*M3*L4*L5+L3*AA4*L2*M4* L5-AA1*M1*w;


    Lsim1[5][0].Cos = -AA4*AA4*M2*L3*M5-D_d2*D_d2*M2*L3*M5+2*D_d3*M5*M2*D_d4+2*D_d3*M5*M2*L4*D_d5+2*AA3*AA2*L2*M3*M5+2*AA3*L5*M2*M3*AA5+D_d4*D_d4*M2*L3*M5+2*L4*D_d5*D_d4*M2* L3*M5+AA2*AA2*M2*L3*M5+D_d3*D_d3*M2*L3*M5-2*AA2*L5*L2*L3*AA5+D_d5*D_d5*M2*L3*M5+AA3*AA3*M2*L3*M5+AA5*AA5*M2*L3*M5;

    Lsim1[5][0].Sin = 2*D_d4*AA2*L2*M5+2*AA2*L3*M5*L2*D_d3+2*L4*D_d5*AA2*L2 *M5+2*D_d5*L5*M2*L4*AA5+2*AA2*L3*M5*D_d2+2*L3*AA5* D_d3*M2*L5-2*D_d3*M3*M5*M2*AA3+2*AA5*D_d4*M2*L5-2*AA4 *M5*M2*M4*D_d5;

    Lsim1[5][0].Const = -AA4*AA4*L2*M3*M5-2*D_d3*M3*M5*D_d2+2*L4*D_d5*D_d4*L2* M3*M5+2*M3*AA5*AA2*M2*L5+D_d4*D_d4*L2*M3* M5+ AA5*AA5*L2*M3*M5+2*AA2*L3*M5*M2*AA3-2*AA3*L5*L2*L3*AA5+AA2*AA2*L2*M3*M5-D_d2*D_d2*L2*M3* M5-D_d3*D_d3*L2*M3*M5+AA3*AA3*L2*M3*M5+D_d5*D_d5*L2*M3*M5;

    Lsim1[5][1].Cos = 2*D_d4*M4*M5*M2*L3*AA4-2*AA5*D_d3*M2*L4*L5+2*AA4*D_d3 *M2*M4*M5-2*AA2*L4*M5*L2*D_d3+2*D_d5*AA3* M2*M3*M5-2*AA2*L4*M5*L2*L3*D_d4-2*AA2*L4*M5*D_d2-2*D_d4* L4*L5*M2*L3*AA5-2*D_d5*AA2*L2*L3*M5 -2*D_d5*L5*M2*L3*AA5+2*AA3*L4*M5*M2*M3*D_d4;
    Lsim1[5][1].Sin =  2*D_d4*M5*M2*D_d5+D_d4*D_d4*M2*L4*M5+2*AA5*AA4*M2*M4* L5+2*D_d3*L3*M5*M2*D_d5+2*AA2*M4*M5*L2*AA4 +AA2*AA2* M2*L4*M5+D_d3*D_d3*M2*L4*M5+D_d5*D_d5*M2*L4*M5+2*L3*D_d4*D_d3*M2*L4*M5+AA5*AA5*M2*L4*M5-D_d2*D_d2*M2*L4*M5-AA3*AA3*M2*L4*M5+AA4*AA4*M2*L4* M5-2*AA2*L4*L5*L2*AA5;
    Lsim1[5][1].Const = -2*AA3*L4*M5*L2*D_d3-2*AA3*L4*M5*L2*L3*D_d4+2*D_d4* M4*M5*L2*M3*AA4+2*M3*D_d4*AA2*M2*L4*M5-2 *D_d5 *L5*L2*M3*AA5-2*D_d5*AA3*L2*L3*M5-2*D_d4*L4*L5*L2*M3*AA5-2*AA3*L4*M5*D_d2+2*AA2*M3*M5*M2*D_d5;

    Lsim1[5][2].Cos =  -2*D_d3*L5*M2*AA5+2*AA3*M5*M2*M3*D_d4+2*AA4*M5*M2* L3*M4*D_d5-2*AA2*M5*L2*L3*D_d4-2*AA5*D_d4*M2*L3*L5+2*AA3*M5*M2*M3*L4*D_d5-2*AA2*M5*L2*L3*L4*D_d5-2*AA2*M5*D_d2-2*AA2*M5*L2*D_d3-2*D_d5*L5*M2*L3*L4*AA5;
    Lsim1[5][2].Sin =  D_d4*D_d4*M2*M5+2*L3*L4*D_d5*D_d3*M2*M5-2*AA5*AA2*L2* L5+AA5*AA5*M2*M5-AA3*AA3*M2*M5+D_d3*D_d3*M2* M5-AA4*AA4 *M2*M5+2*L3*D_d4*D_d3*M2*M5+AA2*AA2*M2*M5-D_d2*D_d2* M2*M5+2*L4*D_d5*D_d4*M2*M5+D_d5*D_d5*M2*M5;
    Lsim1[5][2].Const = 2*AA4*M5*L2*M3*M4*D_d5-2*AA3*M5*L2*D_d3+2*M3*L4*D_d5*AA2*M2*M5-2*AA3*M5*L2*L3*D_d4-2*AA3*M5*D_d2+2*M3*D_d4*AA2*M2*M5-2*D_d5*L5*L2*M3*L4*AA5-2*AA3*M5 *L2*L3*L4*D_d5-2*AA5*D_d4*L2*M3*L5;


    Lsim1[5][3].Cos =  -D_d3*D_d3*M2*L3*L4*M5+2*AA5*AA2*L2*L3*L4*L5- 2*AA5*AA4*M2*L3*M4*L5+D_d2*D_d2*M2*L3*L4*M5-2 *AA3*AA2*L2*M3*L4*M5-2*D_d4*M5*M2*L3*D_d5-D_d4*D_d4* M2*L3*L4*M5-2*D_d3*L4*M5*M2*D_d4-AA3*AA3*M2 * L3*L4*M5-AA4*AA4*M2*L3*L4*M5-2*D_d5*D_d3*M2* M5-2*AA5*AA3*M2*M3*L4*L5-2*AA4*AA2*L2*L3*M4 *M5-D_d5*D_d5*M2*L3*L4*M5+2*AA4*AA3*M2*M3*M4* M5-AA2*AA2*M2*L3*L4*M5-AA5*AA5*M2*L3*L4*M5;
    Lsim1[5][3].Sin = 2*D_d4*M4*M5*M2*AA4-2*AA2*L3*L4*M5*L2*D_d3-2*D_d5* L5*M2*AA5+2*D_d3*L3*M4*M5*M2*AA4-2*D_d3*L3* L4 *L5*M2*AA5-2*AA2*L3*L4*M5*D_d2+2*D_d3*M3*L4*M5*M2*AA3-2*D_d4*AA2*L2*L4*M5-2*AA2*M5*L2*D_d5 -2*D_d4*L4*L5*M2*AA5;
    Lsim1[5][3].Const =   D_d3*D_d3*L2*M3*L4*M5+D_d2*D_d2*L2*M3*L4*M5+2* AA5*AA3*L2*L3*L4*L5-AA3*AA3*L2*M3*L4*M5-2* AA2*M3*L4*L5*M2*AA5-2*AA4*AA3*L2*L3*M4*M5- AA4*AA4*L2*M3*L4*M5-2*AA5*AA4*L2*M3*M4*L5-2 *AA2*L3*L4*M5*M2*AA3-AA5*AA5*L2*M3*L4*M5+2* D_d3*M3*L4*M5*D_d2+2*AA2*M3*M4*M5*M2*AA4-2*D_d4*M5 *L2*M3*D_d5-D_d4*D_d4*L2*M3*L4*M5-D_d5*D_d5*L2*M3*L4*M5-AA2*AA2*L2*M3*L4*M5;

    Lsim1[5][4].Cos = -2*AA2*M4*L5*L2*D_d3-2*D_d4*L4*L5*M2*L3*AA4+2*D_d3* M5*M2*M4*AA5-2*AA2*M4*L5*D_d2+2*AA3*M4*L5 *M2 *M3*D_d4-2*AA4*D_d3*M2*L4*L5-2*AA2*M4*L5*L2*L3*D_d4+2*M4*AA5*D_d4*M2*L3*M5-2*D_d5*L5*M2*L3*AA4;
    Lsim1[5][4].Sin = D_d3*D_d3*M2*M4*L5+2*L3*D_d4*D_d3*M2*M4*L5+AA4*AA4*M2*M4*L5+2*AA4*M5*M2*L4*AA5-D_d2*D_d2*M2*M4*L5+2*M4*AA5*AA2*L2*M5+AA2*AA2*M2*M4*L5-D_d5*D_d5*M2 *M4*L5+AA5*AA5*M2*M4*L5+D_d4*D_d4*M2*M4*L5- 2*AA2*L4*L5*L2*AA4-AA3*AA3*M2*M4*L5;
    Lsim1[5][4].Const = -2*AA3*M4*L5*D_d2-2*D_d5*L5*L2*M3*AA4+2*M3*D_d4*AA2*M2*M4*L5-2*AA3*M4*L5*L2*L3*D_d4-2*AA3*M4*L5*L2*D_d3+2*M4*AA5*D_d4*L2*M3*M5-2*D_d4*L4*L5*L2*M3*AA4;

    Lsim1[5][5].Cos = -D_d3*D_d3*M2*L3*M4*L5-D_d4*D_d4*M2*L3*M4*L5-AA5*AA5*M2*L3*M4*L5+D_d5*D_d5*M2*L3*M4*L5+D_d2*D_d2*M2*L3*M4*L5-2*AA2*M5*L2*L3*M4*AA5-2*AA3*AA2*L2*M3*M4*L5-2*D_d3*M4*L5*M2*D_d4-AA4*AA4*M2*L3*M4*L5-AA2*AA2*M2*L3*M4*L5-AA3*AA3*M2* L3*M4*L5+2*AA4*AA2*L2*L3*L4*L5+2*AA3*M5* M2*M3*M4*AA5-2*AA4*M5*M2*L3*L4*AA5-2*AA4*AA3*M2*M3*L4*L5;
    Lsim1[5][5].Sin = 2*L3*M4*AA5*D_d3*M2*M5-2*D_d3*L3*L4*L5*M2*AA4-2*D_d4*L4*L5*M2*AA4+2*M4*AA5*D_d4*M2*M5-2*AA2*L3*M4*L5*L2*D_d3-2*AA2*L3*M4*L5*D_d2+2*D_d3*M3*M4*L5*M2*AA3-2*D_d5*L5*M2*AA4-2*D_d4*AA2*L2*M4*L5;
    Lsim1[5][5].Const = -AA3*AA3*L2*M3*M4*L5-2*AA2*L3*M4*L5*M2*AA3-2*AA3*M5*L2*L3*M4*AA5-(AA42+AA52)*L2*M3*M4*L5+2*D_d3*M3*M4*L5*D_d2+D_d2*D_d2*L2*M3*M4*L5-2*AA4*M5*L2*M3*L4*AA5+2*AA4*AA3*L2*L3*L4*L5-AA2*AA2*L2*M3*M4*L5+2*M3*M4*AA5*AA2*M2*M5-2*AA2*M3*L4*L5*M2*AA4+D_d3*D_d3*L2*M3*M4*L5+D_d5*D_d5*L2*M3*M4*L5-D_d42*L2*M3*M4*L5;

    Lsim1[5][6].Cos = 2*M4*D_d5*AA2*L2*M3*M5+2*AA4*M5*M2*M3*D_d4+2*D_d5*L5*M2*M3*M4*AA5+2*AA4*M5*M2*M3*L4*D_d5+2* M4*D_d5*AA3*M2*L3*M5;
    Lsim1[5][6].Sin = -2*AA4*M5*M2*AA3-2*D_d3*M3*M5*M2*M4*D_d5;
    Lsim1[5][6].Const = -2*AA4*M5*D_d2-2*AA4*M5*L2*D_d3+2*M4*D_d5*AA3*L2*M3*M5+2*AA2*L3*M5*M2*M4*D_d5-2*AA4*M5*L2*L3*D_d4-2*D_d5* L5*L2*L3*M4*AA5-2*AA4*M5*L2*L3*L4*D_d5;

    Lsim1[5][7].Cos =  AA3*AA3*M2*M3*M4*M5-D_d2*D_d2*M2*M3*M4*M5+D_d3*D_d3*M2*M3*M4*M5-2*AA5*AA4*M2*M3*L4*L5-2*AA3*M4*L5*M2*L3*AA5+D_d5*D_d5*M2*M3*M4*M5-2*AA3*L4*M5*M2*L3*AA4+AA5*AA5*M2*M3*M4*M5-D_d4*D_d4*M2*M3*M4*M5-2*AA3*AA2*L2*L3*M4*M5+AA4*AA4*M2*M3*M4*M5-2*AA2*L4*M5*L2*M3*AA4+AA2*AA2*M2*M3*M4*M5-2*AA2*M4*L5*L2*M3*AA5;
    Lsim1[5][7].Sin = 2*D_d3*L3*M4*M5*M2*AA3+2*AA2*M3*M4*M5*D_d2+2*M3*AA5*D_d3*M2*M4*L5+2*M3*AA4*D_d3*M2*L4*M5+2*AA2*M3*M4*M5*L2*D_d3+2*D_d4*M4*M5*M2*AA3;
    Lsim1[5][7].Const = -D_d5*D_d5*L2*L3*M4*M5-AA4*AA4*L2*L3*M4*M5-2*L3*AA4*AA2*M2*L4*M5+2*D_d3*L3*M4*M5*D_d2+D_d3*D_d3*L2*L3*M4*M5-2*AA3*L4*M5*L2*M3*AA4+D_d2*D_d2*L2*L3*M4*M5-2*AA3*M4*L5*L2*M3*AA5+D_d4*D_d4*L2*L3*M4*M5-AA3*AA3*L2*L3*M4*M5-AA5*AA5*L2*L3*M4*M5+2*D_d4*M4*M5*D_d2-AA2*AA2*L2*L3*M4*M5+2*AA5*AA4*L2*L3*L4*L5+2*AA2*M3*M4*M5*M2*AA3+2*D_d4*M4*M5*L2*D_d3-2*L3*AA5*AA2*M2*M4*L5;

    Lsim1[5][8].Cos = -2*AA3*M4*L5*M2*L3*AA4-2*L4*AA5*AA3*M2*L3*M5-AA5*AA5*M2*M3*L4*L5-2*AA2*M4*L5*L2*M3*AA4 +2*AA3*AA2*L2*L3*L4*L5-AA4*AA4*M2*M3*L4*L5+ D_d2*D_d2*M2*M3*L4*L5+2*D_d5*L5*M2*M3*D_d4-D_d3*D_d3*M2*M3*L4*L5-AA32*M2*M3*L4*L5+2*AA4*M5 *M2*M3*M4*AA5+D_d52*M2*M3*L4*L5-2*L4*AA5* AA2*L2*M3*M5-AA2*AA2*M2*M3*L4*L5+D_d4*D_d4*M2* M3*L4*L5;
    Lsim1[5][8].Sin = -2*D_d3*L3*L4*L5*M2*AA3-2*AA2*M3*L4*L5*D_d2-2*D_d4* L4*L5*M2*AA3+2*M3*AA4*D_d3*M2*M4*L5-2*D_d5*L5 *M2*AA3-2*AA2*M3*L4*L5*L2*D_d3+2*D_d3*M3*M5*M2*L4*AA5;
    Lsim1[5][8].Const = AA5*AA5*L2*L3*L4*L5+AA2*AA2*L2*L3*L4*L5-D_d3*D_d3*L2*L3*L4*L5-D_d4*D_d4*L2*L3*L4*L5-2*D_d3* L3*L4*L5*D_d2-2*D_d4*L4*L5*D_d2+2*p*u*L1*r+r*r*L1 *w-2*p*u*L1*D_d1-2*q*v*L1*D_d1-2*r*D_d1*L1*w-AA1*AA1*L1*w+2*q*v*L1*r+D_d1*D_d1*L1*w-p*p*L1*w-q*q*L1*w-2*D_d5*L5*L2*D_d3-2* D_d5*L5*D_d2-2*D_d4*L4*L5*L2*D_d3-2*L4*AA5*AA3*L2*M3* M5-2*D_d5*L5*L2*L3*D_d4-2*AA2*M3*L4*L5*M2*AA3 -2*AA1*v*M1*p+2*q*AA1*M1*u-2*AA2*L3*M5*M2*L4*AA5-D_d5*D_d5*L2*L3*L4*L5+AA4*AA4*L2*L3*L4*L5-D_d2*D_d2* L2*L3*L4*L5-2*L3*AA4*AA2*M2*M4*L5-2*AA4* M5*L2*L3*M4*AA5-2*AA3*M4*L5*L2*M3*AA4+AA3*AA3*L2*L3*L4*L5;





    Lsim2[0][0].Cos = -L4 * AA5;
    Lsim2[0][0].Sin = 0.0;
    Lsim2[0][0].Const = 0.0;

    Lsim2[0][1].Cos = 0.0;
    Lsim2[0][1].Sin = -L3 * AA5;
    Lsim2[0][1].Const = 0.0;

    Lsim2[0][2].Cos = 0.0;
    Lsim2[0][2].Sin = -L3*L4*AA5;
    Lsim2[0][2].Const = 0.0;

    Lsim2[0][3].Cos = AA5;
    Lsim2[0][3].Sin = 0.0;
    Lsim2[0][3].Const = 0.0;

    Lsim2[0][4].Cos = M4*D_d5;
    Lsim2[0][4].Sin = -L3 * AA4;
    Lsim2[0][4].Const = 0.0;

    Lsim2[0][5].Cos = AA4;
    Lsim2[0][5].Sin = L3*M4*D_d5;
    Lsim2[0][5].Const = 0.0;

    Lsim2[0][6].Cos = 0.0;
    Lsim2[0][6].Sin = M3*M4*AA5;
    Lsim2[0][6].Const = 0.0;

    Lsim2[0][7].Cos = 0.0;
    Lsim2[0][7].Sin = 0.0;
    Lsim2[0][7].Const = 0.0;

    Lsim2[0][8].Cos = AA3;
    Lsim2[0][8].Sin = M3*L4*D_d5+D_d4*M3;
    Lsim2[0][8].Const = AA2;


    Lsim2[1][0].Cos = 0.0;
    Lsim2[1][0].Sin = L2*L4*AA5;
    Lsim2[1][0].Const = 0.0;

    Lsim2[1][1].Cos = -L2*L3*AA5;
    Lsim2[1][1].Sin = 0.0;
    Lsim2[1][1].Const = M2*M3*AA5;

    Lsim2[1][2].Cos = -L2*L3*L4*AA5;
    Lsim2[1][2].Sin = 0.0;
    Lsim2[1][2].Const = M2*M3*L4*AA5;

    Lsim2[1][3].Cos = 0.0;
    Lsim2[1][3].Sin = -L2*AA5;
    Lsim2[1][3].Const = 0.0;

    Lsim2[1][4].Cos = -AA4*L3*L2;
    Lsim2[1][4].Sin = -L2*M4*D_d5;
    Lsim2[1][4].Const = M2*M3*AA4;

    Lsim2[1][5].Cos =  L2*L3*M4*D_d5;
    Lsim2[1][5].Sin = -L2*AA4;
    Lsim2[1][5].Const = -M2*M3*M4*D_d5;

    Lsim2[1][6].Cos = L2*M3*M4*AA5;
    Lsim2[1][6].Sin = 0.0;
    Lsim2[1][6].Const = M2*L3*M4*AA5;

    Lsim2[1][7].Cos = 0.0;
    Lsim2[1][7].Sin = 0.0;
    Lsim2[1][7].Const = 0.0;

    Lsim2[1][8].Cos = L2*M3*D_d4+L2*M3*L4*D_d5;
    Lsim2[1][8].Sin = -L2*AA3;
    Lsim2[1][8].Const = M2*L3*D_d4+M2*D_d3+M2*L3*L4*D_d5;

    Lsim2[2][0].Cos = 0.0;
    Lsim2[2][0].Sin = -L3*M5;
    Lsim2[2][0].Const = 0.0;

    Lsim2[2][1].Cos = L4*M5;
    Lsim2[2][1].Sin = 0.0;
    Lsim2[2][1].Const = 0.0;

    Lsim2[2][2].Cos = M5;
    Lsim2[2][2].Sin = 0.0;
    Lsim2[2][2].Const = 0.0;

    Lsim2[2][3].Cos = 0.0;
    Lsim2[2][3].Sin = L3*L4*M5;
    Lsim2[2][3].Const = 0.0;

    Lsim2[2][4].Cos = M4*L5;
    Lsim2[2][4].Sin = 0.0;
    Lsim2[2][4].Const = 0.0;

    Lsim2[2][5].Cos = 0.0;
    Lsim2[2][5].Sin = L3*M4*L5;
    Lsim2[2][5].Const = 0.0;

    Lsim2[2][6].Cos = 0.0;
    Lsim2[2][6].Sin = 0.0;
    Lsim2[2][6].Const = 0.0;

    Lsim2[2][7].Cos = 0.0;
    Lsim2[2][7].Sin = -M3*M4*M5;
    Lsim2[2][7].Const = 0.0;

    Lsim2[2][8].Cos = 0.0;
    Lsim2[2][8].Sin = M3*L4*L5;
    Lsim2[2][8].Const = 0.0;

    Lsim2[3][0].Cos = -L2*L3*M5;
    Lsim2[3][0].Sin = 0.0;
    Lsim2[3][0].Const = M2*M3*M5;

    Lsim2[3][1].Cos = 0.0;
    Lsim2[3][1].Sin = -L2*L4*M5;
    Lsim2[3][1].Const = 0.0;

    Lsim2[3][2].Cos = 0.0;
    Lsim2[3][2].Sin = -L2*M5;
    Lsim2[3][2].Const = 0.0;

    Lsim2[3][3].Cos = L2*L3*L4*M5;
    Lsim2[3][3].Sin = 0.0;
    Lsim2[3][3].Const = -M2*M3*L4*M5;

    Lsim2[3][4].Cos = 0.0;
    Lsim2[3][4].Sin = -L2*M4*L5;
    Lsim2[3][4].Const = 0.0;

    Lsim2[3][5].Cos = L2*L3*M4*L5;
    Lsim2[3][5].Sin = 0.0;
    Lsim2[3][5].Const = -M2*M3*M4*L5;

    Lsim2[3][6].Cos = 0.0;
    Lsim2[3][6].Sin = 0.0;
    Lsim2[3][6].Const = 0.0;

    Lsim2[3][7].Cos = -L2*M3*M4*M5;
    Lsim2[3][7].Sin = 0.0;
    Lsim2[3][7].Const = -M2*L3*M4*M5;

    Lsim2[3][8].Cos = L2*M3*L4*L5;
    Lsim2[3][8].Sin = 0.0;
    Lsim2[3][8].Const = M2*L3*L4*L5;

    Lsim2[4][0].Cos = D_d4*M5+D_d2*L2*L3*M5+L4*D_d5*M5+D_d3*L3*M5; 
    Lsim2[4][0].Sin = L3*AA5*L5-AA3*M3*M5;
    Lsim2[4][0].Const = -D_d2*M2*M3*M5;

    Lsim2[4][1].Cos = -AA5*L4*L5+AA4*M4*M5;
    Lsim2[4][1].Sin = D_d3*L4*M5+D_d2*L2*L4*M5+D_d5*L3*M5+L3*D_d4*L4*M5;
    Lsim2[4][1].Const = 0.0;

    Lsim2[4][2].Cos = -AA5*L5;
    Lsim2[4][2].Sin = L3*D_d4*M5+D_d2*L2*M5+D_d3*M5+L3*L4*D_d5*M5;
    Lsim2[4][2].Const = 0.0;

    Lsim2[4][3].Cos = -D_d5*M5-D_d2*L2*L3*L4*M5-D_d3*L3*L4*M5-D_d4*L4*M5;
    Lsim2[4][3].Sin = -AA5*L3*L4*L5+AA3*M3*L4*M5+AA4*L3*M4*M5;
    Lsim2[4][3].Const = D_d2*M2*M3*L4*M5;

    Lsim2[4][4].Cos = M4*AA5*M5-AA4*L4*L5;
    Lsim2[4][4].Sin = L3*D_d4*M4*L5+D_d2*L2*M4*L5+D_d3*M4*L5;
    Lsim2[4][4].Const = 0.0;

    Lsim2[4][5].Cos = -D_d2*L2*L3*M4*L5-D_d3*L3*M4*L5-D_d4*M4*L5;
    Lsim2[4][5].Sin = AA3*M3*M4*L5+L3*M4*AA5*M5-AA4*L3*L4*L5;
    Lsim2[4][5].Const =  D_d2*M2*M3*M4*L5;

    Lsim2[4][6].Cos = 0.0;
    Lsim2[4][6].Sin = -M4*D_d5*M3*M5;
    Lsim2[4][6].Const = 0.0;

    Lsim2[4][7].Cos = D_d3*M3*M4*M5+D_d2*L2*M3*M4*M5;
    Lsim2[4][7].Sin = AA3*L3*M4*M5+M3*AA5*M4*L5+M3*AA4*L4*M5;
    Lsim2[4][7].Const = D_d2*M2*L3*M4*M5;

    Lsim2[4][8].Cos = -D_d2*L2*M3*L4*L5-D_d3*M3*L4*L5;
    Lsim2[4][8].Sin = L4*AA5*M3*M5+M3*AA4*M4*L5-AA3*L3*L4*L5;
    Lsim2[4][8].Const = -D_d2*M2*L3*L4*L5;
   

    Lsim2[5][0].Cos = -AA2*M2*L3*M5-AA3*L2*M3*M5+L2*L3*AA5*L5;
    Lsim2[5][0].Sin = -D_d4*L2*M5-L4*D_d5*L2*M5-L2*D_d3*L3*M5-D_d2*L3*M5;
    Lsim2[5][0].Const = -AA2*L2*M3*M5-M2*AA3*L3*M5-M3*AA5*M2*L5;

    Lsim2[5][1].Cos = L2*L3*D_d4*L4*M5+D_d2*L4*M5+L2*D_d3*L4*M5+D_d5*L2*L3*M5;
    Lsim2[5][1].Sin = -L2*AA4*M4*M5+L2*AA5*L4*L5-AA2*M2*L4*M5;
    Lsim2[5][1].Const = -M3*D_d4*M2*L4*M5-M2*D_d5*M3*M5;

    Lsim2[5][2].Cos = L2*L3*L4*D_d5*M5+L2*L3*D_d4*M5+L2*D_d3*M5+D_d2*M5;
    Lsim2[5][2].Sin = -AA2*M2*M5+AA5*L2*L5;
    Lsim2[5][2].Const = -M3*D_d4*M2*M5-M3*L4*D_d5*M2*M5;

    Lsim2[5][3].Cos = -AA5*L2*L3*L4*L5+AA2*M2*L3*L4*M5+AA3*L2*M3*L4*M5+AA4*L2*L3*M4*M5;
    Lsim2[5][3].Sin = D_d2*L3*L4*M5+L2*D_d3*L3*L4*M5+L2*D_d5*M5+D_d4*L2*L4*M5;
    Lsim2[5][3].Const = AA2*L2*M3*L4*M5+M2*AA5*M3*L4*L5-M2*AA4*M3*M4*M5+M2*AA3*L3*L4*M5;

    Lsim2[5][4].Cos =  L2*D_d3*M4*L5+D_d2*M4*L5+L2*L3*D_d4*M4*L5;
    Lsim2[5][4].Sin = -M4*AA5*L2*M5+L2*AA4*L4*L5-AA2*M2*M4*L5;
    Lsim2[5][4].Const = -M3*D_d4*M2*M4*L5;

    Lsim2[5][5].Cos = L2*L3*M4*AA5*M5+AA2*M2*L3*M4*L5-AA4*L2*L3*L4*L5+AA3*L2*M3*M4*L5;
    Lsim2[5][5].Sin = D_d4*L2*M4*L5+D_d2*L3*M4*L5+L2*D_d3*L3*M4*L5;
    Lsim2[5][5].Const = M2*AA4*M3*L4*L5+AA2*L2*M3*M4*L5+M2*AA3*L3*M4*L5-M3*M4*AA5*M2*M5;

    Lsim2[5][6].Cos = -M4*D_d5*L2*M3*M5;
    Lsim2[5][6].Sin = 0.0;
    Lsim2[5][6].Const = -M2*M4*D_d5*L3*M5;

    Lsim2[5][7].Cos = L2*M3*AA5*M4*L5-AA2*M2*M3*M4*M5+L2*M3*AA4*L4*M5+AA3*L2*L3*M4*M5;
    Lsim2[5][7].Sin = -D_d2*M3*M4*M5-L2*D_d3*M3*M4*M5;
    Lsim2[5][7].Const = L3*AA5*M2*M4*L5+L3*AA4*M2*L4*M5+AA2*L2*L3*M4*M5-M2*AA3*M3*M4*M5;

    Lsim2[5][8].Cos = L4*AA5*L2*M3*M5+L2*M3*AA4*M4*L5-AA3*L2*L3*L4*L5+AA2*M2*M3*L4*L5;
    Lsim2[5][8].Sin = D_d2*M3*L4*L5+L2*D_d3*M3*L4*L5;
    Lsim2[5][8].Const = M2*AA3*M3*L4*L5-AA2*L2*L3*L4*L5+L3*AA4*M2*M4*L5+M2*L4*AA5*L3*M5;
   

    Lsim2[6][0].Cos =  -2*AA4*M5*M4*D_d5+2*AA5*D_d4*L5-2*D_d3*M3*M5*AA3+2*D_d5*L5*L4*AA5-2*D_d2*L2*M3*M5*AA3-2*D_d2*M2*L3*M5*AA2+2*L2*L3*AA5*D_d2*L5+2*L3*AA5*D_d3*L5;
    Lsim2[6][0].Sin = -2*D_d2*L2*M5*L4*D_d5-D_d52*L3*M5-D_d42*L3*M5-2*D_d3*M5*D_d4-2*L2*D_d3*D_d2*L3*M5-AA52*L3*M5-2*L4*D_d5*D_d4*L3*M5-2*AA3*L5*M3*AA5-D_d32*L3*M5-2*D_d2*L2*M5*D_d4-D_d22*L3*M5-AA32*L3*M5-2*D_d3*M5*L4*D_d5+AA42*L3*M5+AA22*L3*M5;
    Lsim2[6][0].Const = -2*D_d3*M3*M5*AA2-2*M2*AA3*D_d2*L3*M5-2*D_d2*M2*L5*M3*AA5-2*D_d2*L2*M3*M5*AA2;

    Lsim2[6][1].Cos = AA42*L4*M5+D_d22*L4*M5+2*L2*L3*D_d4*D_d2*L4*M5+AA52*L4*M5+D_d52*L4*M5+2*L2*D_d3*D_d2*L4*M5+D_d42*L4*M5-AA22*L4*M5-AA32*L4*M5+2*D_d3*L3*M5*D_d5+2*L3*D_d4*D_d3*L4*M5+2*D_d2*L2*L3*M5*D_d5+2*D_d4*M5*D_d5+2*AA5*AA4*M4*L5+D_d32*L4*M5;
    Lsim2[6][1].Sin = -2*AA3*L4*M5*M3*D_d4-2*D_d4*M4*M5*L3*AA4+2*L2*AA5*D_d2*L4*L5-2*L2*AA4*D_d2*M4*M5-2*D_d2*M2*L4*M5*AA2+2*D_d4*L4*L5*L3*AA5-2*D_d5*AA3*M3*M5+2*D_d5*L5*L3*AA5+2*AA5*D_d3*L4*L5-2*AA4*D_d3*M4*M5;
    Lsim2[6][1].Const = -2*M2*D_d5*D_d2*M3*M5-2*AA3*L4*M5*AA2-2*D_d2*M2*L4*M5*M3*D_d4;

    Lsim2[6][2].Cos = 2*L4*D_d5*D_d4*M5+2*L3*L4*D_d5*D_d3*M5+2*L2*L3*L4*D_d5*D_d2*M5+2*L2*L3*D_d4*D_d2*M5+2*L2*D_d3*D_d2*M5+2*L3*D_d4*D_d3*M5+D_d22*M5+D_d42*M5+AA52*M5+D_d52*M5-AA22*M5-AA32*M5-AA42*M5+D_d32*M5;
    Lsim2[6][2].Sin = -2*D_d2*M2*M5*AA2-2*AA4*M5*L3*M4*D_d5-2*AA3*M5*M3*D_d4-2*AA3*M5*M3*L4*D_d5+2*D_d2*L2*L5*AA5+2*D_d3*L5*AA5+2*AA5*D_d4*L3*L5+2*D_d5*L5*L3*L4*AA5;
    Lsim2[6][2].Const = -2*D_d2*M2*M5*M3*D_d4-2*D_d2*M2*M5*M3*L4*D_d5-2*AA3*M5*AA2;

    Lsim2[6][3].Cos = -2*D_d4*L4*L5*AA5-2*D_d5*L5*AA5-2*D_d3*L3*L4*L5*AA5+2*D_d2*L2*M3*L4*M5*AA3+2*D_d2*L2*L3*M4*M5* AA4+2*D_d2*M2*L3*L4*M5*AA2+2*D_d3*M3*L4*M5*AA3+2*D_d4*M4*M5*AA4-2*D_d2*L2*L3*L4*L5*AA5+2*D_d3*L3*M4*M5*AA4;
    Lsim2[6][3].Sin = 2*AA5*AA4*L3*M4*L5-2*AA4*AA3*M3*M4*M5+2*D_d2*L2*L4*M5*D_d4+2*AA5*AA3*M3*L4*L5+2*D_d3*L4*M5*D_d4+2*L2*D_d3*D_d2*L3*L4*M5+AA52*L3*L4*M5+D_d32*L3*L4*M5+AA32*L3*L4*M5+2*L2*D_d5*D_d2*M5 +D_d52*L3*L4*M5+D_d42*L3*L4*M5+AA42*L3*L4*M5+2*D_d4*M5*L3*D_d5+D_d22*L3*L4*M5-AA22*L3*L4*M5+2*D_d5*D_d3*M5;
    Lsim2[6][3].Const = 2*M2*AA5*D_d2*M3*L4*L5+2*D_d3*M3*L4*M5*AA2-2*M2*AA4*D_d2*M3*M4*M5+2*D_d2*L2*M3*L4*M5*AA2+2*M2*AA3*D_d2*L3*L4*M5;

    Lsim2[6][4].Cos = AA52*M4*L5+D_d32*M4*L5+2*AA4*M5*L4*AA5-D_d52*M4*L5+AA42*M4*L5-AA22*M4*L5+D_d42*M4*L5+D_d22*M4*L5+2*L2*D_d3*D_d2*M4*L5+2*L2*L3*D_d4*D_d2*M4*L5+2*L3*D_d4*D_d3*M4*L5-AA32*M4*L5;
    Lsim2[6][4].Sin = -2*D_d2*L2*M5*M4*AA5-2*D_d3*M5*M4*AA5-2*AA3*M4*L5*M3*D_d4-2*D_d2*M2*M4*L5*AA2+2*L2*AA4*D_d2*L4*L5-2*M4*AA5*D_d4*L3*M5+2*D_d4*L4*L5*L3*AA4+2*AA4*D_d3*L4*L5+2*D_d5*L5*L3*AA4;
    Lsim2[6][4].Const = -2*AA3*M4*L5*AA2-2*D_d2*M2*M4*L5*M3*D_d4;

    Lsim2[6][5].Cos = 2*D_d3*M3*M4*L5*AA3+2*L3*M4*AA5*D_d3*M5+2*L2*L3*M4*AA5*D_d2*M5-2*D_d2*L2*L3*L4*L5*AA4-2*D_d3*L3*L4*L5*AA4+2*M4*AA5*D_d4*M5+2*D_d2*M2*L3*M4*L5*AA2-2*D_d4*L4*L5*AA4-2*D_d5*L5*AA4+2*D_d2*L2*M3*M4*L5*AA3;
    Lsim2[6][5].Sin = -2*AA3*M5*M3*M4*AA5-AA22*L3*M4*L5+D_d32*L3*M4*L5+2*D_d2*L2*M4*L5*D_d4+2*D_d3*M4*L5*D_d4+2*AA4*M5*L3*L4*AA5+AA32*L3*M4*L5+AA42*L3*M4*L5+AA52*L3*M4*L5-D_d52*L3*M4*L5+D_d22*L3*M4*L5+2*AA4*AA3*M3*L4*L5+2*L2*D_d3*D_d2*L3*M4*L5+D_d42*L3*M4*L5;
    Lsim2[6][5].Const = -2*D_d2*M2*M5*M3*M4*AA5+2*M2*AA3*D_d2*L3*M4*L5+2*D_d2*L2*M3*M4*L5*AA2+2*M2*AA4*D_d2*M3*L4*L5+2*D_d3*M3*M4*L5*AA2;

    Lsim2[6][6].Cos = -2*D_d2*L2*M3*M5*M4*D_d5-2*D_d3*M3*M5*M4*D_d5-2*AA4*M5*AA3;
    Lsim2[6][6].Sin = -2*AA4*M5*M3*D_d4-2*AA4*M5*M3*L4*D_d5-2*M4*D_d5*AA3*L3*M5-2*D_d5*L5*M3*M4*AA5;
    Lsim2[6][6].Const = -2*AA4*M5*AA2-2*M2*M4*D_d5*D_d2*L3*M5;

    Lsim2[6][7].Cos = 2*D_d2*L2*L3*M4*M5*AA3+2*D_d4*M4*M5*AA3+2*M3*AA5*D_d3*M4*L5+2*D_d3*L3*M4*M5*AA3+2*L2*M3*AA5*D_d2*M4*L5-2*D_d2*M2*M3*M4*M5*AA2+2*L2*M3*AA4*D_d2*L4*M5+2*M3*AA4*D_d3*L4*M5;
    Lsim2[6][7].Sin =  -2*L2*D_d3*D_d2*M3*M4*M5-D_d22*M3*M4*M5-AA52*M3*M4*M5-D_d52*M3*M4*M5-D_d32*M3*M4*M5+D_d42*M3*M4*M5-AA42*M3*M4*M5+2*AA3*M4*L5*L3*AA5+2*AA3*L4*M5*L3*AA4+2*AA5*AA4*M3*L4*L5+AA22*M3*M4*M5-AA32*M3*M4*M5;
    Lsim2[6][7].Const = -2*M2*AA3*D_d2*M3*M4*M5+2*D_d2*M2*L4*M5*L3*AA4+2*D_d3*L3*M4*M5*AA2+2*D_d2*M2*M4*L5*L3*AA5+2*D_d4*M4*M5*AA2+2*D_d2*L2*L3*M4*M5*AA2;

    Lsim2[6][8].Cos = 2*D_d3*M3*M5*L4*AA5-2*D_d5*L5*AA3+2*D_d2*L2*M3*M5*L4*AA5+2*D_d2*M2*M3*L4*L5*AA2-2*D_d4*L4*L5*AA3-2*D_d2*L2*L3*L4*L5*AA3+2*L2*M3*AA4*D_d2*M4*L5+2*M3*AA4*D_d3*M4*L5-2*D_d3*L3*L4*L5*AA3;
    Lsim2[6][8].Sin = -2*D_d5*L5*M3*D_d4-AA22*M3*L4*L5+2*L2*D_d3*D_d2*M3*L4*L5-D_d52*M3*L4*L5-D_d42*M3*L4*L5+ AA32*M3*L4*L5+2*AA3*M4*L5*L3*AA4+2*L4*AA5*AA3*L3*M5+AA52*M3*L4*L5-2*AA4*M5*M3*M4*AA5 +D_d32*M3*L4*L5+AA42*M3*L4*L5+D_d22*M3*L4*L5;
    Lsim2[6][8].Const = -2*D_d2*L2*L3*L4*L5*AA2-2*D_d5*L5*AA2-2*D_d4*L4*L5*AA2+2*D_d2*M2*M4*L5*L3*AA4-2*D_d3*L3*L4*L5* AA2+2*M2*AA3*D_d2*M3*L4*L5+2*M2*L4*AA5*D_d2*L3*M5;
   

    Lsim2[7][0].Cos = -2*D_d3*M5*L2*D_d4-2*D_d3*M5*L2*L4*D_d5-2*L4*D_d5*D_d4*L2*L3*M5-2*AA2*L5*M2*L3*AA5-AA32*L2*L3*M5-AA52*L2*L3*M5-D_d32*L2*L3*M5-D_d42*L2*L3*M5-AA22*L2*L3*M5-2*AA3*L5*L2*M3*AA5-D_d52*L2*L3*M5+AA42*L2*L3*M5-2*D_d2*L3*M5*D_d3-D_d22*L2*L3*M5-2*D_d4*D_d2*M5-2*L4*D_d5*D_d2*M5+2*AA3*AA2*M2*M3*M5;
    Lsim2[7][0].Sin = 2*AA4*M5*L2*M4*D_d5+2*D_d4*AA2*M2*M5-2*AA5*D_d4*L2*L5+2*D_d3*M3*M5*L2*AA3+2*L4*D_d5*AA2*M2*M5+2*AA2*L3*M5*M2*D_d3+2*AA3*D_d2*M3*M5-2*L3*AA5*D_d3*L2*L5-2*D_d5*L5*L2*L4*AA5-2*D_d2*L5*L3*AA5;
    Lsim2[7][0].Const = D_d52*M2*M3*M5-2*AA2*L3*M5*L2*AA3+D_d22*M2*M3*M5+AA22*M2*M3*M5-2*AA3*L5*M2*L3*AA5-AA42*M2*M3*M5+AA52*M2*M3*M5-2*M3*AA5*AA2*L2*L5-D_d32*M2*M3*M5+D_d42*M2*M3*M5+2*L4*D_d5*D_d4*M2*M3*M5+AA32*M2*M3*M5;

    Lsim2[7][1].Cos = -2*AA3*L4*M5*L2*M3*D_d4-2*D_d2*M4*M5*AA4-2*AA4*D_d3*L2*M4*M5-2*D_d5*AA3*L2*M3*M5+2*AA5*D_d3*L2*L4*L5+2*D_d2*L4*L5*AA5+2*D_d4*L4*L5*L2*L3*AA5-2*D_d5*AA2*M2*L3*M5-2*AA2*L4*M5*M2*L3*D_d4-2*AA2*L4*M5*M2*D_d3+2*D_d5*L5*L2*L3*AA5-2*D_d4*M4*M5*L2*L3*AA4;
    Lsim2[7][1].Sin = AA32*L2*L4*M5-2*AA2*L4*L5*M2*AA5-AA22*L2*L4*M5-2*D_d3*L3*M5*L2*D_d5-2*D_d2*L4*M5*L3*D_d4-D_d22*L2*L4*M5-2*D_d5*D_d2*L3*M5-2*D_d4*M5*L2*D_d5-2*AA5*AA4*L2*M4*L5-D_d42*L2*L4*M5-AA52*L2*L4*M5+2*AA2*M4*M5*M2*AA4-D_d32*L2*L4*M5-2*L3*D_d4*D_d3*L2*L4*M5-AA42*L2*L4*M5-D_d52*L2*L4*M5-2*D_d2*L4*M5*D_d3;
    Lsim2[7][1].Const = -2*AA3*L4*M5*M2*D_d3+2*D_d4*M4*M5*M2*M3*AA4-2*M3*D_d4*AA2*L2*L4*M5-2*D_d4*L4*L5*M2*M3*AA5-2*D_d5*AA3*M2*L3*M5-2*AA2*M3*M5*L2*D_d5-2*D_d5*L5*M2*M3*AA5-2*AA3*L4*M5*M2*L3*D_d4;

    Lsim2[7][2].Cos = 2*D_d5*L5*L2*L3*L4*AA5+2*AA5*D_d2*L5-2*AA2*M5*M2*L3*D_d4-2*AA2*M5*M2*D_d3-2*AA3*M5*L2*M3*L4*D_d5+2*D_d3*L5*L2*AA5-2*AA2*M5*M2*L3*L4*D_d5-2*AA4*M5*L2*L3*M4*D_d5+2*AA5*D_d4*L2*L3*L5-2*AA3*M5*L2*M3*D_d4;
    Lsim2[7][2].Sin = -AA52*L2*M5-D_d42*L2*M5-2*D_d2*M5*D_d3-2*D_d2*M5*L3*L4*D_d5-2*L3*L4*D_d5*D_d3*L2*M5-2*L3*D_d4*D_d3*L2*M5-AA22*L2*M5-D_d52*L2*M5-D_d32*L2*M5+AA32*L2*M5-2*AA5*AA2*M2*L5-2*D_d2*M5*L3*D_d4-2*L4*D_d5*D_d4*L2*M5+AA42*L2*M5-D_d22*L2*M5;
    Lsim2[7][2].Const = -2*AA5*D_d4*M2*M3*L5-2*AA3*M5*M2*D_d3-2*D_d5*L5*M2*M3*L4*AA5-2*AA3*M5*M2*L3*D_d4+2*AA4*M5*M2*M3*M4*D_d5-2*M3*L4*D_d5*AA2*L2*M5-2*AA3*M5*M2*L3*L4*D_d5-2*M3*D_d4*AA2*L2*M5;

    Lsim2[7][3].Cos = 2*AA5*AA3*L2*M3*L4*L5+2*AA5*AA2*M2*L3*L4*L5+2*D_d5*D_d3*L2*M5-2*AA4*AA3*L2*M3*M4*M5+2*D_d4*D_d2*L4*M5-2*AA4*AA2*M2*L3*M4*M5+2*D_d2*L3*L4*M5*D_d3+2*D_d2*M5*D_d5+AA52*L2*L3*L4*M5+D_d32*L2*L3*L4*M5+AA32*L2*L3*L4*M5+D_d52*L2*L3*L4*M5+D_d42*L2*L3*L4*M5+AA22*L2*L3*L4*M5+2*AA5*AA4*L2*L3*M4*L5+2*D_d3*L4*M5*L2*D_d4+AA42*L2*L3*L4*M5+2*D_d4*M5*L2*L3*D_d5+D_d22*L2*L3*L4*M5-2*AA3*AA2*M2*M3*L4*M5;
    Lsim2[7][3].Sin = 2*D_d4*L4*L5*L2*AA5-2*AA2*L3*L4*M5*M2*D_d3-2*AA4*D_d2*L3*M4*M5-2*D_d3*M3*L4*M5*L2*AA3-2*AA3*D_d2*M3*L4*M5+2*D_d3*L3*L4*L5*L2*AA5+2*AA5*D_d2*L3*L4*L5-2*AA2*M5*M2*D_d5-2*D_d4*AA2*M2*L4*M5-2*D_d4*M4*M5*L2*AA4+2*D_d5*L5*L2*AA5-2*D_d3*L3*M4*M5*L2*AA4;
    Lsim2[7][3].Const = -2*AA4*AA3*M2*L3*M4*M5-AA42*M2*M3*L4*M5+2*AA2*M3*L4*L5*L2*AA5-D_d22*M2*M3*L4*M5-D_d52*M2*M3*L4*M5-AA32*M2*M3*L4*M5+2*AA5*AA3*M2*L3*L4*L5-2*AA5*AA4*M2*M3*M4*L5-2*AA2*M3*M4*M5*L2*AA4-D_d42*M2*M3*L4*M5-AA52*M2*M3*L4*M5-AA22*M2*M3*L4*M5+2*AA2*L3*L4*M5*L2*AA3-2*D_d4*M5*M2*M3*D_d5+D_d32*M2*M3*L4*M5;

    Lsim2[7][4].Cos = -2*M4*AA5*D_d4*L2*L3*M5-2*AA2*M4*L5*M2*D_d3-2*M4*AA5*D_d2*M5+2*D_d5*L5*L2*L3*AA4-2*AA2*M4*L5*M2*L3*D_d4-2*D_d3*M5*L2*M4*AA5+2*D_d4*L4*L5*L2*L3*AA4+2*D_d2*L4*L5*AA4-2*AA3*M4*L5*L2*M3*D_d4+2*AA4*D_d3*L2*L4*L5;
    Lsim2[7][4].Sin = -2*D_d2*M4*L5*D_d3-AA22*L2*M4*L5-AA42*L2*M4*L5-AA52*L2*M4*L5+2*M4*AA5*AA2*M2*M5-2*D_d2*M4*L5*L3*D_d4-2*L3*D_d4*D_d3*L2*M4*L5-2*AA4*M5*L2*L4*AA5-D_d42*L2*M4*L5+D_d52*L2*M4*L5+AA32*L2*M4*L5-D_d22*L2*M4*L5-D_d32*L2*M4*L5-2*AA2*L4*L5*M2*AA4;
    Lsim2[7][4].Const = 2*M4*AA5*D_d4*M2*M3*M5-2*M3*D_d4*AA2*L2*M4*L5-2*AA3*M4*L5*M2*L3*D_d4-2*AA3*M4*L5*M2*D_d3-2*D_d4*L4*L5*M2*M3*AA4-2*D_d5*L5*M2*M3*AA4;

    Lsim2[7][5].Cos = AA22*L2*L3*M4*L5-D_d52*L2*L3*M4*L5+2*D_d4*D_d2*M4*L5+2*AA4*M5*L2*L3*L4*AA5-2*AA3*AA2*M2*M3*M4*L5+D_d32*L2*L3*M4*L5+AA52*L2*L3*M4*L5+2*D_d2*L3*M4*L5*D_d3+2*AA4*AA2*M2*L3*L4*L5+2*D_d3*M4*L5*L2*D_d4-2*AA2*M5*M2*L3*M4*AA5+2*AA4*AA3*L2*M3*L4*L5-2*AA3*M5*L2*M3*M4*AA5+AA32*L2*L3*M4*L5+D_d42*L2*L3*M4*L5+D_d22*L2*L3*M4*L5+AA42*L2*L3*M4*L5;
    Lsim2[7][5].Sin =  2*D_d4*L4*L5*L2*AA4+2*D_d5*L5*L2*AA4-2*D_d2*M5*L3*M4*AA5-2*L3*M4*AA5*D_d3*L2*M5+2*AA4*D_d2*L3*L4*L5+2*D_d3*L3*L4*L5*L2*AA4-2*AA2*L3*M4*L5*M2*D_d3-2*D_d4*AA2*M2*M4*L5-2*AA3*D_d2*M3*M4*L5-2*D_d3*M3*M4*L5*L2*AA3-2*M4*AA5*D_d4*L2*M5;
    Lsim2[7][5].Const = D_d52*M2*M3*M4*L5+2*AA2*L3*M4*L5*L2*AA3-2*M3*M4*AA5*AA2*L2*M5-D_d42*M2*M3*M4*L5-AA52*M2*M3*M4*L5-AA22*M2*M3*M4*L5+2*AA4*AA3*M2*L3*L4*L5-AA32*M2*M3*M4*L5-AA42*M2*M3*M4*L5-D_d22*M2*M3*M4*L5-2*AA3*M5*M2*L3*M4*AA5+D_d32*M2*M3*M4*L5+2*AA2*M3*L4*L5*L2*AA4-2*AA4*M5*M2*M3*L4*AA5;

    Lsim2[7][6].Cos = -2*M4*D_d5*AA3*L2*L3*M5-2*AA4*M5*L2*M3*L4*D_d5-2*D_d5*L5*L2*M3*M4*AA5+2*M4*D_d5*AA2*M2*M3*M5-2*AA4*M5*L2*M3*D_d4;
    Lsim2[7][6].Sin = 2*M4*D_d5*D_d2*M3*M5+2*D_d3*M3*M5*L2*M4*D_d5+2*AA4*M5*L2*AA3;
    Lsim2[7][6].Const = -2*AA4*M5*M2*L3*D_d4-2*AA2*L3*M5*L2*M4*D_d5-2*D_d5*L5*M2*L3*M4*AA5+2*M4*D_d5*AA3*M2*M3*M5-2*AA4*M5*M2*D_d3-2*AA4*M5*M2*L3*L4*D_d5;

    Lsim2[7][7].Cos = -2*AA2*M4*L5*M2*M3*AA5+D_d42*L2*M3*M4*M5-D_d32*L2*M3*M4*M5-D_d22*L2*M3*M4*M5-2*AA2*L4*M5*M2*M3*AA4-AA42*L2*M3*M4*M5-AA22*L2*M3*M4*M5-2*AA3*AA2*M2*L3*M4*M5-2*D_d2*M3*M4*M5*D_d3+2*AA5*AA4*L2*M3*L4*L5-D_d52*L2*M3*M4*M5+2*AA3*L4*M5*L2*L3*AA4+2*AA3*M4*L5*L2*L3*AA5-AA32*L2*M3*M4*M5-AA52*L2*M3*M4*M5;
    Lsim2[7][7].Sin = -2*AA3*D_d2*L3*M4*M5-2*D_d3*L3*M4*M5*L2*AA3-2*D_d4*M4*M5*L2*AA3-2*D_d2*M4*L5*M3*AA5-2*M3*AA4*D_d3*L2*L4*M5+2*AA2*M3*M4*M5*M2*D_d3-2*D_d2*L4*M5*M3*AA4-2*M3*AA5*D_d3*L2*M4*L5;
    Lsim2[7][7].Const = D_d42*M2*L3*M4*M5+D_d32*M2*L3*M4*M5-D_d52*M2*L3*M4*M5-2*AA2*M3*M4*M5*L2*AA3+2*D_d4*M4*M5*M2*D_d3-AA52*M2*L3*M4*M5-AA42*M2*L3*M4*M5-AA22*M2*L3*M4*M5-2*AA3*L4*M5*M2*M3*AA4+2*AA5*AA4*M2*L3*L4*L5+2*L3*AA5*AA2*L2*M4*L5-D_d22*M2*L3*M4*M5-AA32*M2*L3*M4*M5-2*AA3*M4*L5*M2*M3*AA5+2*L3*AA4*AA2*L2*L4*M5;

    Lsim2[7][8].Cos = -D_d42*L2*M3*L4*L5-2*AA4*M5*L2*M3*M4*AA5+2*AA3*AA2*M2*L3*L4*L5+2*AA3*M4*L5*L2*L3*AA4-2*D_d5*L5*L2*M3*D_d4+2*L4*AA5*AA3*L2*L3*M5-D_d52*L2*M3*L4*L5+D_d22*L2*M3*L4*L5+AA42*L2*M3*L4*L5+2*D_d2*M3*L4*L5*D_d3+AA32*L2*M3*L4*L5+AA52*L2*M3*L4*L5+AA22*L2*M3*L4*L5+D_d32*L2*M3*L4*L5-2*AA2*M4*L5*M2*M3*AA4-2*L4*AA5*AA2*M2*M3*M5;
    Lsim2[7][8].Sin = -2*L4*AA5*D_d2*M3*M5+2*D_d3*L3*L4*L5*L2*AA3-2*D_d2*M4*L5*M3*AA4+2*D_d5*L5*L2*AA3-2*D_d3*M3*M5*L2*L4*AA5-2*AA2*M3*L4*L5*M2*D_d3+2*D_d4*L4*L5*L2*AA3-2*M3*AA4*D_d3*L2*M4*L5+2*AA3*D_d2*L3*L4*L5;
    Lsim2[7][8].Const = -D_d52*M2*L3*L4*L5-2*AA4*M5*M2*L3*M4*AA5+2*AA2*L3*M5*L2*L4*AA5+2*AA2*M3*L4*L5*L2*AA3+AA52*M2*L3*L4*L5-2*AA3*M4*L5*M2*M3*AA4-2*D_d5*L5*M2*L3*D_d4+2*L3*AA4*AA2*L2*M4*L5-D_d32*M2*L3*L4*L5+AA32*M2*L3*L4*L5-D_d42*M2*L3*L4*L5-2*D_d5*L5*M2*D_d3+AA22*M2*L3*L4*L5-2*D_d4*L4*L5*M2*D_d3-2*L4*AA5*AA3*M2*M3*M5+AA42*M2*L3*L4*L5+D_d22*M2*L3*L4*L5;

}
