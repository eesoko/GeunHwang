/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * minOrMax.c
 *
 * Code generation for function 'minOrMax'
 *
 */

/* Include files */
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"

/* Function Definitions */
double maximum(const double x[6], int *idx)
{
  double ex;
  int b_idx;
  int b_k;
  if (!rtIsNaN(x[0])) {
    b_idx = 1;
  } else {
    int k;
    boolean_T exitg1;
    b_idx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 7)) {
      if (!rtIsNaN(x[k - 1])) {
        b_idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }
  if (b_idx == 0) {
    ex = x[0];
    *idx = 1;
  } else {
    ex = x[b_idx - 1];
    *idx = b_idx;
    b_idx++;
    for (b_k = b_idx; b_k < 7; b_k++) {
      double d;
      d = x[b_k - 1];
      if (ex < d) {
        ex = d;
        *idx = b_k;
      }
    }
  }
  return ex;
}

/* End of code generation (minOrMax.c) */
