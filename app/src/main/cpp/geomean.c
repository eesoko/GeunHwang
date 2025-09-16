/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * geomean.c
 *
 * Code generation for function 'geomean'
 *
 */

/* Include files */
#include "geomean.h"
#include "blockedSummation.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "omp.h"
#include <math.h>

/* Function Definitions */
double geomean(const emxArray_real_T *x)
{
  emxArray_real_T *b_x;
  const double *x_data;
  double d;
  double m;
  double *b_x_data;
  int b_loop_ub;
  int k;
  int loop_ub;
  x_data = x->data;
  emxInit_real_T(&b_x, 1);
  loop_ub = x->size[0];
  b_loop_ub = b_x->size[0];
  b_x->size[0] = x->size[0];
  emxEnsureCapacity_real_T(b_x, b_loop_ub);
  b_x_data = b_x->data;
  b_loop_ub = x->size[0];
  if (x->size[0] < 1200) {
    for (k = 0; k < loop_ub; k++) {
      b_x_data[k] = log(x_data[k]);
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(d)

    for (k = 0; k < b_loop_ub; k++) {
      d = log(x_data[k]);
      b_x_data[k] = d;
    }
  }
  m = exp(colMajorFlatIter(b_x, b_x->size[0]) / (double)x->size[0]);
  emxFree_real_T(&b_x);
  return m;
}

/* End of code generation (geomean.c) */
