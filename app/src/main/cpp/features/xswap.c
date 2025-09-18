/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xswap.c
 *
 * Code generation for function 'xswap'
 *
 */

/* Include files */
#include "xswap.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_xswap(int n, double x_data[], int ix0, int iy0)
{
  int i;
  int k;
  i = (unsigned char)n;
  for (k = 0; k < i; k++) {
    double temp;
    int i1;
    int temp_tmp;
    temp_tmp = (ix0 + k) - 1;
    temp = x_data[temp_tmp];
    i1 = (iy0 + k) - 1;
    x_data[temp_tmp] = x_data[i1];
    x_data[i1] = temp;
  }
}

void xswap(double x[9], int ix0, int iy0)
{
  double temp;
  temp = x[ix0 - 1];
  x[ix0 - 1] = x[iy0 - 1];
  x[iy0 - 1] = temp;
  temp = x[ix0];
  x[ix0] = x[iy0];
  x[iy0] = temp;
  temp = x[ix0 + 1];
  x[ix0 + 1] = x[iy0 + 1];
  x[iy0 + 1] = temp;
}

/* End of code generation (xswap.c) */
