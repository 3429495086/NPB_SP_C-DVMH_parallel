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
// step in the z-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the z-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void z_solve()
{
  int i, j, k, k1, k2, m;
  double ru1, fac1, fac2;

  if (timeron) timer_start(t_zsolve);
  #pragma dvm region
  {
    #pragma dvm parallel(2) private(k, k1, k2, m, ru1, fac1, fac2)
    for (j = 1; j <= ny2; j++) {
      for (i = 1; i <= nx2; i++) {
        //---------------------------------------------------------------------
        //location part
        //---------------------------------------------------------------------
        double lhs_loc[PROBLEM_SIZE][5];
        double lhsp_loc[PROBLEM_SIZE][5];
        double lhsm_loc[PROBLEM_SIZE][5];
        double cv[PROBLEM_SIZE];
        double rhos[PROBLEM_SIZE];

        int nj = nz2+1;
        for (m = 0; m < 5; m++) {
          lhs_loc [0][m] = 0.0;
          lhsp_loc[0][m] = 0.0;
          lhsm_loc[0][m] = 0.0;
          lhs_loc [nj][m] = 0.0;
          lhsp_loc[nj][m] = 0.0;
          lhsm_loc[nj][m] = 0.0;
        }
        lhs_loc [0][2] = 1.0;
        lhsp_loc[0][2] = 1.0;
        lhsm_loc[0][2] = 1.0;
        lhs_loc [nj][2] = 1.0;
        lhsp_loc[nj][2] = 1.0;
        lhsm_loc[nj][2] = 1.0;
        //---------------------------------------------------------------------
        // Computes the left hand side for the three z-factors   
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
        // first fill the lhs for the u-eigenvalue                          
        //---------------------------------------------------------------------
        for (k = 0; k <= nz2+1; k++) {
          ru1 = c3c4*rho_i[k][j][i];
          cv[k] = ws[k][j][i];
          rhos[k] = max(max(dz4+con43*ru1, dz5+c1c5*ru1), max(dzmax+ru1, dz1));
        }

        for (k = 1; k <= nz2; k++) {
          lhs_loc[k][0] =  0.0;
          lhs_loc[k][1] = -dttz2 * cv[k-1] - dttz1 * rhos[k-1];
          lhs_loc[k][2] =  1.0 + c2dttz1 * rhos[k];
          lhs_loc[k][3] =  dttz2 * cv[k+1] - dttz1 * rhos[k+1];
          lhs_loc[k][4] =  0.0;
        }
        //---------------------------------------------------------------------
        // add fourth order dissipation                                  
        //---------------------------------------------------------------------

        k = 1;
        lhs_loc[k][2] += comz5;
        lhs_loc[k][3] -= comz4;
        lhs_loc[k][4] += comz1;
        lhs_loc[k+1][1] -= comz4;
        lhs_loc[k+1][2] += comz6;
        lhs_loc[k+1][3] -= comz4;
        lhs_loc[k+1][4] += comz1;

        for (k = 3; k <= grid_points[2]-4; k++) {
          lhs_loc[k][0] += comz1;
          lhs_loc[k][1] -= comz4;
          lhs_loc[k][2] += comz6;
          lhs_loc[k][3] -= comz4;
          lhs_loc[k][4] += comz1;
        }

        k = nz2-1;
        lhs_loc[k][0] += comz1;
        lhs_loc[k][1] -= comz4;
        lhs_loc[k][2] += comz6;
        lhs_loc[k][3] -= comz4;
        lhs_loc[k+1][0] += comz1;
        lhs_loc[k+1][1] -= comz4;
        lhs_loc[k+1][2] += comz5;
        //---------------------------------------------------------------------
        // subsequently, fill the other factors (u+c), (u-c) 
        //---------------------------------------------------------------------
        for (k = 1; k <= grid_points[2]-2; k++) {
          lhsp_loc[k][0] = lhs_loc[k][0];
          lhsp_loc[k][1] = lhs_loc[k][1] - dttz2 * speed[k-1][j][i];
          lhsp_loc[k][2] = lhs_loc[k][2];
          lhsp_loc[k][3] = lhs_loc[k][3] + dttz2 * speed[k+1][j][i];
          lhsp_loc[k][4] = lhs_loc[k][4];
          lhsm_loc[k][0] = lhs_loc[k][0];
          lhsm_loc[k][1] = lhs_loc[k][1] + dttz2 * speed[k-1][j][i];
          lhsm_loc[k][2] = lhs_loc[k][2];
          lhsm_loc[k][3] = lhs_loc[k][3] - dttz2 * speed[k+1][j][i];
          lhsm_loc[k][4] = lhs_loc[k][4];
        }
        //---------------------------------------------------------------------
        // FORWARD ELIMINATION  
        //---------------------------------------------------------------------
        for (k = 0; k <= grid_points[2]-3; k++) {
          k1 = k + 1;
          k2 = k + 2;

          // LHS
          fac1 = 1.0/lhs_loc[k][2];
          lhs_loc[k][3] *= fac1;
          lhs_loc[k][4] *= fac1;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] *= fac1;
          }
          lhs_loc[k1][2] -= lhs_loc[k1][1]*lhs_loc[k][3];
          lhs_loc[k1][3] -= lhs_loc[k1][1]*lhs_loc[k][4];
          for (m = 0; m < 3; m++) {
            rhs[k1][j][i][m] -= lhs_loc[k1][1]*rhs[k][j][i][m];
          }
          lhs_loc[k2][1] -= lhs_loc[k2][0]*lhs_loc[k][3];
          lhs_loc[k2][2] -= lhs_loc[k2][0]*lhs_loc[k][4];
          for (m = 0; m < 3; m++) {
            rhs[k2][j][i][m] -= lhs_loc[k2][0]*rhs[k][j][i][m];
          }
          //---------------------------------------------------------------------
          // for (the u+c and the u-c factors               
          //---------------------------------------------------------------------
          //m = 3;
          fac1 = 1.0/lhsp_loc[k][2];
          lhsp_loc[k][3]   *= fac1;
          lhsp_loc[k][4]   *= fac1;
          rhs[k][j][i][3]  *= fac1;
          lhsp_loc[k1][2]  -= lhsp_loc[k1][1]*lhsp_loc[k][3];
          lhsp_loc[k1][3]  -= lhsp_loc[k1][1]*lhsp_loc[k][4];
          rhs[k1][j][i][3] -= lhsp_loc[k1][1]*rhs[k][j][i][3];
          lhsp_loc[k2][1]  -= lhsp_loc[k2][0]*lhsp_loc[k][3];
          lhsp_loc[k2][2]  -= lhsp_loc[k2][0]*lhsp_loc[k][4];
          rhs[k2][j][i][3] -= lhsp_loc[k2][0]*rhs[k][j][i][3];

          //m = 4;
          fac1 = 1.0/lhsm_loc[k][2];
          lhsm_loc[k][3]   *= fac1;
          lhsm_loc[k][4]   *= fac1;
          rhs[k][j][i][4]  *= fac1;
          lhsm_loc[k1][2]  -= lhsm_loc[k1][1]*lhsm_loc[k][3];
          lhsm_loc[k1][3]  -= lhsm_loc[k1][1]*lhsm_loc[k][4];
          rhs[k1][j][i][4] -= lhsm_loc[k1][1]*rhs[k][j][i][4];
          lhsm_loc[k2][1]  -= lhsm_loc[k2][0]*lhsm_loc[k][3];
          lhsm_loc[k2][2]  -= lhsm_loc[k2][0]*lhsm_loc[k][4];
          rhs[k2][j][i][4] -= lhsm_loc[k2][0]*rhs[k][j][i][4];
        }
        //---------------------------------------------------------------------
        // The last two rows in this grid block are a bit different, 
        // since they for (not have two more rows available for the
        // elimination of off-diagonal entries
        //---------------------------------------------------------------------
        k  = grid_points[2]-2;
        k1 = grid_points[2]-1;
        fac1 = 1.0/lhs_loc[k][2];
        lhs_loc[k][3] *= fac1;
        lhs_loc[k][4] *= fac1;
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs_loc[k1][2] -= lhs_loc[k1][1]*lhs_loc[k][3];
        lhs_loc[k1][3] -= lhs_loc[k1][1]*lhs_loc[k][4];
        for (m = 0; m < 3; m++) {
          rhs[k1][j][i][m] -= lhs_loc[k1][1]*rhs[k][j][i][m];
        }

        //---------------------------------------------------------------------
        // scale the last row immediately
        //---------------------------------------------------------------------
        fac2 = 1.0/lhs_loc[k1][2];
        for (m = 0; m < 3; m++) {
          rhs[k1][j][i][m] = fac2*rhs[k1][j][i][m];
        }
        //---------------------------------------------------------------------
        // And again the last two rows separately
        //---------------------------------------------------------------------
        //m = 3;
        fac1 = 1.0/lhsp_loc[k][2];
        lhsp_loc[k][3]   *= fac1;
        lhsp_loc[k][4]   *= fac1;
        rhs[k][j][i][3]  *= fac1;
        lhsp_loc[k1][2]  -= lhsp_loc[k1][1]*lhsp_loc[k][3];
        lhsp_loc[k1][3]  -= lhsp_loc[k1][1]*lhsp_loc[k][4];
        rhs[k1][j][i][3] -= lhsp_loc[k1][1]*rhs[k][j][i][3];
        rhs[k1][j][i][3] /= lhsp_loc[k1][2];

        //m = 4;
        fac1 = 1.0/lhsm_loc[k][2];
        lhsm_loc[k][3]   *= fac1;
        lhsm_loc[k][4]   *= fac1;
        rhs[k][j][i][4]  *= fac1;
        lhsm_loc[k1][2]  -= lhsm_loc[k1][1]*lhsm_loc[k][3];
        lhsm_loc[k1][3]  -= lhsm_loc[k1][1]*lhsm_loc[k][4];
        rhs[k1][j][i][4] -= lhsm_loc[k1][1]*rhs[k][j][i][4];
        rhs[k1][j][i][4] /= lhsm_loc[k1][2];
        //---------------------------------------------------------------------
        // BACKSUBSTITUTION 
        //---------------------------------------------------------------------
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] -= lhs_loc[k][3]*rhs[k1][j][i][m];
        }

        rhs[k][j][i][3] -= lhsp_loc[k][3]*rhs[k1][j][i][3];
        rhs[k][j][i][4] -= lhsm_loc[k][3]*rhs[k1][j][i][4];
        //---------------------------------------------------------------------
        // Whether or not this is the last processor, we always have
        // to complete the back-substitution 
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
        // The first three factors
        //---------------------------------------------------------------------
        for (k = grid_points[2]-3; k >= 0; k--) {
          k1 = k + 1;
          k2 = k + 2;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] -= lhs_loc[k][3]*rhs[k1][j][i][m] +
                               lhs_loc[k][4]*rhs[k2][j][i][m];
          }

          //-------------------------------------------------------------------
          // And the remaining two
          //-------------------------------------------------------------------
          rhs[k][j][i][3] -= lhsp_loc[k][3]*rhs[k1][j][i][3] +
                             lhsp_loc[k][4]*rhs[k2][j][i][3];
          rhs[k][j][i][4] -= lhsm_loc[k][3]*rhs[k1][j][i][4] +
                             lhsm_loc[k][4]*rhs[k2][j][i][4];
        }
      }
    }
  }
  if (timeron) timer_stop(t_zsolve);

  tzetar();
}
