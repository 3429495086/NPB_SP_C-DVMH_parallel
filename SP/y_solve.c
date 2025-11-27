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
// step in the y-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the y-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void y_solve()
{
  int i, j, k, j1, j2, m;
  double ru1, fac1, fac2;

  if (timeron) timer_start(t_ysolve);
  #pragma dvm region
  {
    #pragma dvm parallel(2) private(j, j1, j2, m, ru1, fac1, fac2)
    for (k = 1; k <= grid_points[2]-2; k++) {
      for (i = 1; i <= grid_points[0]-2; i++) {

        //---------------------------------------------------------------------
        //location part
        //---------------------------------------------------------------------
        double lhs_loc[PROBLEM_SIZE][5];
        double lhsp_loc[PROBLEM_SIZE][5];
        double lhsm_loc[PROBLEM_SIZE][5];
        double cv[PROBLEM_SIZE];
        double rhoq[PROBLEM_SIZE];

        int nj = ny2+1;
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
        // Computes the left hand side for the three y-factors   
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
        // first fill the lhs for the u-eigenvalue         
        //---------------------------------------------------------------------
        for (j = 0; j <= grid_points[1]-1; j++) {
          ru1 = c3c4*rho_i[k][j][i];
          cv[j] = vs[k][j][i];
          rhoq[j] = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));
        }

        for (j = 1; j <= grid_points[1]-2; j++) {
          lhs_loc[j][0] =  0.0;
          lhs_loc[j][1] = -dtty2 * cv[j-1] - dtty1 * rhoq[j-1];
          lhs_loc[j][2] =  1.0 + c2dtty1 * rhoq[j];
          lhs_loc[j][3] =  dtty2 * cv[j+1] - dtty1 * rhoq[j+1];
          lhs_loc[j][4] =  0.0;
        }
        //---------------------------------------------------------------------
        // add fourth order dissipation                             
        //---------------------------------------------------------------------
        j = 1;
        lhs_loc[j][2] += comz5;
        lhs_loc[j][3] -= comz4;
        lhs_loc[j][4] += comz1;
        lhs_loc[j+1][1] -= comz4;
        lhs_loc[j+1][2] += comz6;
        lhs_loc[j+1][3] -= comz4;
        lhs_loc[j+1][4] += comz1;
        
        for (j = 3; j <= grid_points[1]-4; j++) {
          lhs_loc[j][0] += comz1;
          lhs_loc[j][1] -= comz4;
          lhs_loc[j][2] += comz6;
          lhs_loc[j][3] -= comz4;
          lhs_loc[j][4] += comz1;
        }

        j = grid_points[1]-3;
        lhs_loc[j][0] += comz1;
        lhs_loc[j][1] -= comz4;
        lhs_loc[j][2] += comz6;
        lhs_loc[j][3] -= comz4;
        lhs_loc[j+1][0] += comz1;
        lhs_loc[j+1][1] -= comz4;
        lhs_loc[j+1][2] += comz5;

        for (j = 1; j <= grid_points[1]-2; j++) {
          lhsp_loc[j][0] = lhs_loc[j][0];
          lhsp_loc[j][1] = lhs_loc[j][1] - dtty2 * speed[k][j-1][i];
          lhsp_loc[j][2] = lhs_loc[j][2];
          lhsp_loc[j][3] = lhs_loc[j][3] + dtty2 * speed[k][j+1][i];
          lhsp_loc[j][4] = lhs_loc[j][4];
          lhsm_loc[j][0] = lhs_loc[j][0];
          lhsm_loc[j][1] = lhs_loc[j][1] + dtty2 * speed[k][j-1][i];
          lhsm_loc[j][2] = lhs_loc[j][2];
          lhsm_loc[j][3] = lhs_loc[j][3] - dtty2 * speed[k][j+1][i];
          lhsm_loc[j][4] = lhs_loc[j][4];
        }
        //---------------------------------------------------------------------
        // FORWARD ELIMINATION  
        //---------------------------------------------------------------------
        for (j = 0; j <= grid_points[1]-3; j++) {
          j1 = j + 1;
          j2 = j + 2;
          fac1 = 1.0/lhs_loc[j][2];
          lhs_loc[j][3] *= fac1;
          lhs_loc[j][4] *= fac1;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
          }
          lhs_loc[j1][2] -= lhs_loc[j1][1]*lhs_loc[j][3];
          lhs_loc[j1][3] -= lhs_loc[j1][1]*lhs_loc[j][4];
          for (m = 0; m < 3; m++) {
            rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs_loc[j1][1]*rhs[k][j][i][m];
          }
          lhs_loc[j2][1] -= lhs_loc[j2][0]*lhs_loc[j][3];
          lhs_loc[j2][2] -= lhs_loc[j2][0]*lhs_loc[j][4];
          for (m = 0; m < 3; m++) {
            rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhs_loc[j2][0]*rhs[k][j][i][m];
          }
          //---------------------------------------------------------------------
          // for (the u+c and the u-c factors                 
          //---------------------------------------------------------------------
          //m = 3;
          fac1 = 1.0/lhsp_loc[j][2];
          lhsp_loc[j][3]   *= fac1;
          lhsp_loc[j][4]   *= fac1;
          rhs[k][j][i][3]  *= fac1;
          lhsp_loc[j1][2]  -= lhsp_loc[j1][1]*lhsp_loc[j][3];
          lhsp_loc[j1][3]  -= lhsp_loc[j1][1]*lhsp_loc[j][4];
          rhs[k][j1][i][3] -= lhsp_loc[j1][1]*rhs[k][j][i][3];
          lhsp_loc[j2][1]  -= lhsp_loc[j2][0]*lhsp_loc[j][3];
          lhsp_loc[j2][2]  -= lhsp_loc[j2][0]*lhsp_loc[j][4];
          rhs[k][j2][i][3] -= lhsp_loc[j2][0]*rhs[k][j][i][3];

          //m = 4;
          fac1 = 1.0/lhsm_loc[j][2];
          lhsm_loc[j][3]   *= fac1;
          lhsm_loc[j][4]   *= fac1;
          rhs[k][j][i][4]  *= fac1;
          lhsm_loc[j1][2]  -= lhsm_loc[j1][1]*lhsm_loc[j][3];
          lhsm_loc[j1][3]  -= lhsm_loc[j1][1]*lhsm_loc[j][4];
          rhs[k][j1][i][4] -= lhsm_loc[j1][1]*rhs[k][j][i][4];
          lhsm_loc[j2][1]  -= lhsm_loc[j2][0]*lhsm_loc[j][3];
          lhsm_loc[j2][2]  -= lhsm_loc[j2][0]*lhsm_loc[j][4];
          rhs[k][j2][i][4] -= lhsm_loc[j2][0]*rhs[k][j][i][4];
        }
        //---------------------------------------------------------------------
        // The last two rows in this grid block are a bit different, 
        // since they for (not have two more rows available for the
        // elimination of off-diagonal entries
        //---------------------------------------------------------------------
        j  = grid_points[1]-2;
        j1 = grid_points[1]-1;
        fac1 = 1.0/lhs_loc[j][2];
        lhs_loc[j][3] *= fac1;
        lhs_loc[j][4] *= fac1;
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs_loc[j1][2] -= lhs_loc[j1][1]*lhs_loc[j][3];
        lhs_loc[j1][3] -= lhs_loc[j1][1]*lhs_loc[j][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs_loc[j1][1]*rhs[k][j][i][m];
        }
        //---------------------------------------------------------------------
        // scale the last row immediately 
        //---------------------------------------------------------------------
        fac2 = 1.0/lhs_loc[j1][2];
        for (m = 0; m < 3; m++) {
          rhs[k][j1][i][m] = fac2*rhs[k][j1][i][m];
        }
        //---------------------------------------------------------------------
        // And again the last two rows separately
        //---------------------------------------------------------------------
        //m = 3;
        fac1 = 1.0/lhsp_loc[j][2];
        lhsp_loc[j][3]   *= fac1;
        lhsp_loc[j][4]   *= fac1;
        rhs[k][j][i][3]  *= fac1;
        lhsp_loc[j1][2]  -= lhsp_loc[j1][1]*lhsp_loc[j][3];
        lhsp_loc[j1][3]  -= lhsp_loc[j1][1]*lhsp_loc[j][4];
        rhs[k][j1][i][3] -= lhsp_loc[j1][1]*rhs[k][j][i][3];
        rhs[k][j1][i][3] /= lhsp_loc[j1][2];

        //m = 4;
        fac1 = 1.0/lhsm_loc[j][2];
        lhsm_loc[j][3]   *= fac1;
        lhsm_loc[j][4]   *= fac1;
        rhs[k][j][i][4]  *= fac1;
        lhsm_loc[j1][2]  -= lhsm_loc[j1][1]*lhsm_loc[j][3];
        lhsm_loc[j1][3]  -= lhsm_loc[j1][1]*lhsm_loc[j][4];
        rhs[k][j1][i][4] -= lhsm_loc[j1][1]*rhs[k][j][i][4];
        rhs[k][j1][i][4] /= lhsm_loc[j1][2];

        //---------------------------------------------------------------------
        // BACKSUBSTITUTION 
        //---------------------------------------------------------------------

        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] -= lhs_loc[j][3]*rhs[k][j1][i][m];
        }

        rhs[k][j][i][3] -= lhsp_loc[j][3]*rhs[k][j1][i][3];
        rhs[k][j][i][4] -= lhsm_loc[j][3]*rhs[k][j1][i][4];

        //---------------------------------------------------------------------
        // The first three factors
        //---------------------------------------------------------------------
        for (j = grid_points[1]-3; j >= 0; j--) {
          j1 = j + 1;
          j2 = j + 2;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] -= lhs_loc[j][3]*rhs[k][j1][i][m] +
                               lhs_loc[j][4]*rhs[k][j2][i][m];
          }

          //-------------------------------------------------------------------
          // And the remaining two
          //-------------------------------------------------------------------
          rhs[k][j][i][3] -= lhsp_loc[j][3]*rhs[k][j1][i][3] +
                             lhsp_loc[j][4]*rhs[k][j2][i][3];
          rhs[k][j][i][4] -= lhsm_loc[j][3]*rhs[k][j1][i][4] +
                             lhsm_loc[j][4]*rhs[k][j2][i][4];
        }
      }
    }
  }
     
  if (timeron) timer_stop(t_ysolve);

  pinvr();
}

