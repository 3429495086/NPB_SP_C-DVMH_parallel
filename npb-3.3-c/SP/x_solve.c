//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB SP code. This C        //
//  version is developed by the Center for Manycore Programming at Seoul   //
//  National University and derived from the serial Fortran versions in    //
//  "NPB3.3-SER" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr  //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,    //
//          and Jaejin Lee                                                 //
//-------------------------------------------------------------------------//

#include "header.h"

//---------------------------------------------------------------------
// this function performs the solution of the approximate factorization
// step in the x-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the x-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void x_solve()
{
  int i, j, k, i1, i2, m;
  double ru1, fac1, fac2;

  if (timeron) timer_start(t_xsolve);

  #pragma dvm region 
  {
    #pragma dvm parallel(2) private(i,m,i1,i2,ru1,fac1,fac2)
    for (k = 1; k <= nz2; k++) {
      for (j = 1; j <= ny2; j++) {

        //---------------------------------------------------------------------
        //location part
        //---------------------------------------------------------------------
        double lhs_loc[PROBLEM_SIZE][5];
        double lhsp_loc[PROBLEM_SIZE][5];
        double lhsm_loc[PROBLEM_SIZE][5];
        double cv[PROBLEM_SIZE];
        double rhon[PROBLEM_SIZE];

        int ni = nx2 + 1; 
        for (m = 0; m < 5; m++) {
             lhs_loc[0][m]  = 0.0;
             lhsp_loc[0][m] = 0.0;
             lhsm_loc[0][m] = 0.0;
             lhs_loc[ni][m]  = 0.0;
             lhsp_loc[ni][m] = 0.0;
             lhsm_loc[ni][m] = 0.0;
        }
        lhs_loc[0][2]  = 1.0;
        lhsp_loc[0][2] = 1.0;
        lhsm_loc[0][2] = 1.0;
        lhs_loc[ni][2]  = 1.0;
        lhsp_loc[ni][2] = 1.0;
        lhsm_loc[ni][2] = 1.0;
        //---------------------------------------------------------------------
        // Computes the left hand side for the three x-factors  
        //---------------------------------------------------------------------

        //---------------------------------------------------------------------
        // first fill the lhs for the u-eigenvalue                   
        //---------------------------------------------------------------------
        for (i = 0; i <= grid_points[0]-1; i++) {
          ru1 = c3c4*rho_i[k][j][i];
          cv[i] = us[k][j][i];
          rhon[i] = max(max(dx2+con43*ru1,dx5+c1c5*ru1), max(dxmax+ru1,dx1));
        }

        for (i = 1; i <= nx2; i++) {
          lhs_loc[i][0] =  0.0;
          lhs_loc[i][1] = -dttx2 * cv[i-1] - dttx1 * rhon[i-1];
          lhs_loc[i][2] =  1.0 + c2dttx1 * rhon[i];
          lhs_loc[i][3] =  dttx2 * cv[i+1] - dttx1 * rhon[i+1];
          lhs_loc[i][4] =  0.0;
        }
        //---------------------------------------------------------------------
        // add fourth order dissipation                             
        //---------------------------------------------------------------------
        i = 1;
        lhs_loc[i][2] += comz5;
        lhs_loc[i][3] -= comz4;
        lhs_loc[i][4] += comz1;
        lhs_loc[i+1][1] -= comz4;
        lhs_loc[i+1][2] += comz6;
        lhs_loc[i+1][3] -= comz4;
        lhs_loc[i+1][4] += comz1;

        for (i = 3; i <= grid_points[0]-4; i++) {
          lhs_loc[i][0] += comz1;
          lhs_loc[i][1] -= comz4;
          lhs_loc[i][2] += comz6;
          lhs_loc[i][3] -= comz4;
          lhs_loc[i][4] += comz1;
        }

        i = grid_points[0]-3;
        lhs_loc[i][0] += comz1;
        lhs_loc[i][1] -= comz4;
        lhs_loc[i][2] += comz6;
        lhs_loc[i][3] -= comz4;
        lhs_loc[i+1][0] += comz1;
        lhs_loc[i+1][1] -= comz4;
        lhs_loc[i+1][2] += comz5;
        //---------------------------------------------------------------------
        // subsequently, fill the other factors (u+c), (u-c) by adding to 
        // the first  
        //---------------------------------------------------------------------
        for (i = 1; i <= nx2; i++) {
          lhsp_loc[i][0] = lhs_loc[i][0];
          lhsp_loc[i][1] = lhs_loc[i][1] - dttx2 * speed[k][j][i-1];
          lhsp_loc[i][2] = lhs_loc[i][2];
          lhsp_loc[i][3] = lhs_loc[i][3] + dttx2 * speed[k][j][i+1];
          lhsp_loc[i][4] = lhs_loc[i][4];
          lhsm_loc[i][0] = lhs_loc[i][0];
          lhsm_loc[i][1] = lhs_loc[i][1] + dttx2 * speed[k][j][i-1];
          lhsm_loc[i][2] = lhs_loc[i][2];
          lhsm_loc[i][3] = lhs_loc[i][3] - dttx2 * speed[k][j][i+1];
          lhsm_loc[i][4] = lhs_loc[i][4];
        }
        //---------------------------------------------------------------------
        // FORWARD ELIMINATION  
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
        // perform the Thomas algorithm; first, FORWARD ELIMINATION     
        //---------------------------------------------------------------------
        for (i = 0; i <= grid_points[0]-3; i++) {
          i1 = i + 1;
          i2 = i + 2;
          fac1 = 1.0/lhs_loc[i][2];
          lhs_loc[i][3] *= fac1;
          lhs_loc[i][4] *= fac1;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
          }
          lhs_loc[i1][2] -= lhs_loc[i1][1]*lhs_loc[i][3];
          lhs_loc[i1][3] -= lhs_loc[i1][1]*lhs_loc[i][4];
          for (m = 0; m < 3; m++) {
            rhs[k][j][i1][m] -= lhs_loc[i1][1]*rhs[k][j][i][m];
          }
          lhs_loc[i2][1] -= lhs_loc[i2][0]*lhs_loc[i][3];
          lhs_loc[i2][2] -= lhs_loc[i2][0]*lhs_loc[i][4];
          for (m = 0; m < 3; m++) {
            rhs[k][j][i2][m] -= lhs_loc[i2][0]*rhs[k][j][i][m];
          }
          //---------------------------------------------------------------------
          // for (the u+c and the u-c factors                 
          //---------------------------------------------------------------------
          //m = 3;
          fac1 = 1.0/lhsp_loc[i][2];
          lhsp_loc[i][3]   *= fac1;
          lhsp_loc[i][4]   *= fac1;
          rhs[k][j][i][3]  *= fac1;
          lhsp_loc[i1][2]  -= lhsp_loc[i1][1]*lhsp_loc[i][3];
          lhsp_loc[i1][3]  -= lhsp_loc[i1][1]*lhsp_loc[i][4];
          rhs[k][j][i1][3] -= lhsp_loc[i1][1]*rhs[k][j][i][3];
          lhsp_loc[i2][1]  -= lhsp_loc[i2][0]*lhsp_loc[i][3];
          lhsp_loc[i2][2]  -= lhsp_loc[i2][0]*lhsp_loc[i][4];
          rhs[k][j][i2][3] -= lhsp_loc[i2][0]*rhs[k][j][i][3];

          //m = 4;
          fac1 = 1.0/lhsm_loc[i][2];
          lhsm_loc[i][3]   *= fac1;
          lhsm_loc[i][4]   *= fac1;
          rhs[k][j][i][4]  *= fac1;
          lhsm_loc[i1][2]  -= lhsm_loc[i1][1]*lhsm_loc[i][3];
          lhsm_loc[i1][3]  -= lhsm_loc[i1][1]*lhsm_loc[i][4];
          rhs[k][j][i1][4] -= lhsm_loc[i1][1]*rhs[k][j][i][4];
          lhsm_loc[i2][1]  -= lhsm_loc[i2][0]*lhsm_loc[i][3];
          lhsm_loc[i2][2]  -= lhsm_loc[i2][0]*lhsm_loc[i][4];
          rhs[k][j][i2][4] -= lhsm_loc[i2][0]*rhs[k][j][i][4];
        }
        //---------------------------------------------------------------------
        // And again the last two rows separately
        //---------------------------------------------------------------------
        i  = grid_points[0]-2;
        i1 = grid_points[0]-1;
        fac1 = 1.0/lhs_loc[i][2];
        lhs_loc[i][3] *= fac1;
        lhs_loc[i][4] *= fac1;
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs_loc[i1][2] -= lhs_loc[i1][1]*lhs_loc[i][3];
        lhs_loc[i1][3] -= lhs_loc[i1][1]*lhs_loc[i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i1][m] -= lhs_loc[i1][1]*rhs[k][j][i][m];
        }

        //---------------------------------------------------------------------
        // scale the last row immediately 
        //---------------------------------------------------------------------
        fac2 = 1.0/lhs_loc[i1][2];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i1][m] = fac2*rhs[k][j][i1][m];
        }
        //---------------------------------------------------------------------
        // And again the last two rows separately
        //---------------------------------------------------------------------
        //m = 3;
        fac1 = 1.0/lhsp_loc[i][2];
        lhsp_loc[i][3]   *= fac1;
        lhsp_loc[i][4]   *= fac1;
        rhs[k][j][i][3]  *= fac1;
        lhsp_loc[i1][2]  -= lhsp_loc[i1][1]*lhsp_loc[i][3];
        lhsp_loc[i1][3]  -= lhsp_loc[i1][1]*lhsp_loc[i][4];
        rhs[k][j][i1][3] -= lhsp_loc[i1][1]*rhs[k][j][i][3];
        rhs[k][j][i1][3] /= lhsp_loc[i1][2];

        //m = 4;
        fac1 = 1.0/lhsm_loc[i][2];
        lhsm_loc[i][3]   *= fac1;
        lhsm_loc[i][4]   *= fac1;
        rhs[k][j][i][4]  *= fac1;
        lhsm_loc[i1][2]  -= lhsm_loc[i1][1]*lhsm_loc[i][3];
        lhsm_loc[i1][3]  -= lhsm_loc[i1][1]*lhsm_loc[i][4];
        rhs[k][j][i1][4] -= lhsm_loc[i1][1]*rhs[k][j][i][4];
        rhs[k][j][i1][4] /= lhsm_loc[i1][2];
        //---------------------------------------------------------------------
        // BACKSUBSTITUTION 
        //---------------------------------------------------------------------
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] -= lhs_loc[i][3]*rhs[k][j][i1][m];
        }
        rhs[k][j][i][3] -= lhsp_loc[i][3]*rhs[k][j][i1][3];
        rhs[k][j][i][4] -= lhsm_loc[i][3]*rhs[k][j][i1][4];
        //---------------------------------------------------------------------
        // The first three factors
        //---------------------------------------------------------------------
        for (i = grid_points[0]-3; i >= 0; i--) {
          i1 = i + 1;
          i2 = i + 2;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] -= lhs_loc[i][3]*rhs[k][j][i1][m] +
                               lhs_loc[i][4]*rhs[k][j][i2][m];
          }
          //-------------------------------------------------------------------
          // And the remaining two
          //-------------------------------------------------------------------
          rhs[k][j][i][3] -= lhsp_loc[i][3]*rhs[k][j][i1][3] +
                             lhsp_loc[i][4]*rhs[k][j][i2][3];
          rhs[k][j][i][4] -= lhsm_loc[i][3]*rhs[k][j][i1][4] +
                             lhsm_loc[i][4]*rhs[k][j][i2][4];
        }
      }
    }
  } 
  if (timeron) timer_stop(t_xsolve);
  //---------------------------------------------------------------------
  // Do the block-diagonal inversion          
  //---------------------------------------------------------------------
  ninvr();
}

