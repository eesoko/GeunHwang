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
double b_mean(const emxArray_boolean_T *x)
{
  double y;
  int b_k;
  int k;
  const boolean_T *x_data;
  x_data = x->data;
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    int firstBlockLength;
    int lastBlockLength;
    int nblocks;
    if (x->size[0] <= 1024) {
      firstBlockLength = x->size[0];
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = (int)((unsigned int)x->size[0] >> 10);
      lastBlockLength = x->size[0] - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    y = x_data[0];
    for (k = 2; k <= firstBlockLength; k++) {
      y += (double)x_data[k - 1];
    }
    for (k = 2; k <= nblocks; k++) {
      double bsum;
      int hi;
      firstBlockLength = (k - 1) << 10;
      bsum = x_data[firstBlockLength];
      if (k == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      for (b_k = 2; b_k <= hi; b_k++) {
        bsum += (double)x_data[(firstBlockLength + b_k) - 1];
      }
      y += bsum;
    }
  }
  y /= (double)x->size[0];
  return y;
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
