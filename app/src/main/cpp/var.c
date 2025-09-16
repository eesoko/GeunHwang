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
double var(const double x_data[], int x_size)
{
  emxArray_real_T b_x_data;
  double y;
  int k;
  if (x_size == 0) {
    y = rtNaN;
  } else if (x_size == 1) {
    if ((!rtIsInf(x_data[0])) && (!rtIsNaN(x_data[0]))) {
      y = 0.0;
    } else {
      y = rtNaN;
    }
  } else {
    double xbar;
    b_x_data.data = (double *)&x_data[0];
    b_x_data.size = &x_size;
    b_x_data.allocatedSize = -1;
    b_x_data.numDimensions = 1;
    b_x_data.canFreeData = false;
    xbar = colMajorFlatIter(&b_x_data, x_size) / (double)x_size;
    y = 0.0;
    for (k = 0; k < x_size; k++) {
      double t;
      t = x_data[k] - xbar;
      y += t * t;
    }
    y /= (double)x_size - 1.0;
  }
  return y;
}

/* End of code generation (var.c) */
