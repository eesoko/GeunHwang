/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mean.c
 *
 * Code generation for function 'mean'
 *
 */

/* Include files */
#include "mean.h"
#include "blockedSummation.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
double b_mean(const boolean_T x_data[], int x_size)
{
  int accumulatedData;
  int k;
  if (x_size == 0) {
    accumulatedData = 0;
  } else {
    accumulatedData = x_data[0];
    for (k = 2; k <= x_size; k++) {
      accumulatedData += x_data[k - 1];
    }
  }
  return (double)accumulatedData / (double)x_size;
}

double mean(const emxArray_real_T *x)
{
  double y;
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = colMajorFlatIter(x, x->size[0]);
  }
  y /= (double)x->size[0];
  return y;
}

/* End of code generation (mean.c) */
