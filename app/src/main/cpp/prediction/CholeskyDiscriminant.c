/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CholeskyDiscriminant.c
 *
 * Code generation for function 'CholeskyDiscriminant'
 *
 */

/* Include files */
#include "CholeskyDiscriminant.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void CholeskyDiscriminant(double mah_data[], const int mah_size[2],
                          const boolean_T classLogicalIndices[6],
                          const double X_data[], const int X_size[2],
                          const double Mu[96], const double invD[16],
                          const double invR[256])
{
  double A[16];
  int i;
  int ind;
  int jj;
  int k;
  ind = -1;
  for (jj = 0; jj < 6; jj++) {
    if (classLogicalIndices[jj]) {
      int b_k;
      boolean_T b[16];
      boolean_T exitg1;
      boolean_T y;
      ind++;
      for (k = 0; k < 16; k++) {
        b[k] = rtIsNaN(X_data[X_size[0] * k]);
      }
      y = false;
      b_k = 0;
      exitg1 = false;
      while ((!exitg1) && (b_k < 16)) {
        if (b[b_k]) {
          y = true;
          exitg1 = true;
        } else {
          b_k++;
        }
      }
      if (!y) {
        double X[16];
        double b_y;
        for (k = 0; k < 16; k++) {
          X[k] = (X_data[X_size[0] * k] - Mu[ind + 6 * k]) * invD[k];
        }
        memset(&A[0], 0, sizeof(double) << 4);
        for (k = 0; k < 16; k++) {
          b_y = A[k];
          for (i = 0; i < 16; i++) {
            b_y += X[i] * invR[i + (k << 4)];
          }
          A[k] = b_y * b_y;
        }
        b_y = A[0];
        for (k = 0; k < 15; k++) {
          b_y += A[k + 1];
        }
        mah_data[mah_size[0] * ind] = b_y;
      }
    }
  }
}

/* End of code generation (CholeskyDiscriminant.c) */
