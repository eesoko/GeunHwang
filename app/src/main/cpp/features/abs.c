/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * abs.c
 *
 * Code generation for function 'abs'
 *
 */

/* Include files */
#include "abs.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void b_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  const double *x_data;
  double *y_data;
  int i;
  int k;
  int nx;
  x_data = x->data;
  nx = x->size[0];
  i = y->size[0];
  y->size[0] = x->size[0];
  emxEnsureCapacity_real_T(y, i);
  y_data = y->data;
  for (k = 0; k < nx; k++) {
    y_data[k] = fabs(x_data[k]);
  }
}

/* End of code generation (abs.c) */
