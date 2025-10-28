/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * var.c
 *
 * Code generation for function 'var'
 *
 */

/* Include files */
#include "var.h"
#include "blockedSummation.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"

/* Function Definitions */
double var(const emxArray_real_T *x)
{
  const double *x_data;
  double y;
  int k;
  int n;
  x_data = x->data;
  n = x->size[0];
  if (x->size[0] == 0) {
    y = rtNaN;
  } else if (x->size[0] == 1) {
    if ((!rtIsInf(x_data[0])) && (!rtIsNaN(x_data[0]))) {
      y = 0.0;
    } else {
      y = rtNaN;
    }
  } else {
    double xbar;
    xbar = colMajorFlatIter(x, x->size[0]) / (double)x->size[0];
    y = 0.0;
    for (k = 0; k < n; k++) {
      double t;
      t = x_data[k] - xbar;
      y += t * t;
    }
    y /= (double)x->size[0] - 1.0;
  }
  return y;
}

/* End of code generation (var.c) */
