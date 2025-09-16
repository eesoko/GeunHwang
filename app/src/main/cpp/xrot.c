/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xrot.c
 *
 * Code generation for function 'xrot'
 *
 */

/* Include files */
#include "xrot.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_xrot(int n, double x_data[], int ix0, int iy0, double c, double s)
{
  int k;
  if (n >= 1) {
    int i;
    i = (unsigned char)n;
    for (k = 0; k < i; k++) {
      double b_temp_tmp;
      double temp_tmp;
      int b_temp_tmp_tmp;
      int temp_tmp_tmp;
      temp_tmp_tmp = (iy0 + k) - 1;
      temp_tmp = x_data[temp_tmp_tmp];
      b_temp_tmp_tmp = (ix0 + k) - 1;
      b_temp_tmp = x_data[b_temp_tmp_tmp];
      x_data[temp_tmp_tmp] = c * temp_tmp - s * b_temp_tmp;
      x_data[b_temp_tmp_tmp] = c * b_temp_tmp + s * temp_tmp;
    }
  }
}

void xrot(double x[9], int ix0, int iy0, double c, double s)
{
  double temp;
  double temp_tmp;
  temp = x[iy0 - 1];
  temp_tmp = x[ix0 - 1];
  x[iy0 - 1] = c * temp - s * temp_tmp;
  x[ix0 - 1] = c * temp_tmp + s * temp;
  temp = c * x[ix0] + s * x[iy0];
  x[iy0] = c * x[iy0] - s * x[ix0];
  x[ix0] = temp;
  temp = x[iy0 + 1];
  temp_tmp = x[ix0 + 1];
  x[iy0 + 1] = c * temp - s * temp_tmp;
  x[ix0 + 1] = c * temp_tmp + s * temp;
}

/* End of code generation (xrot.c) */
