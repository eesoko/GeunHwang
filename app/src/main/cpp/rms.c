/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rms.c
 *
 * Code generation for function 'rms'
 *
 */

/* Include files */
#include "rms.h"
#include "blockedSummation.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>

/* Function Definitions */
double rms(const double xIn_data[], int xIn_size)
{
  emxArray_real_T *x;
  double y;
  double *x_data;
  int i;
  int scalarLB;
  int vectorUB;
  emxInit_real_T(&x, 1);
  scalarLB = x->size[0];
  x->size[0] = xIn_size;
  emxEnsureCapacity_real_T(x, scalarLB);
  x_data = x->data;
  scalarLB = (xIn_size / 2) << 1;
  vectorUB = scalarLB - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    __m128d r;
    r = _mm_loadu_pd(&xIn_data[i]);
    _mm_storeu_pd(&x_data[i], _mm_mul_pd(r, r));
  }
  for (i = scalarLB; i < xIn_size; i++) {
    y = xIn_data[i];
    x_data[i] = y * y;
  }
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = colMajorFlatIter(x, x->size[0]);
  }
  y = sqrt(y / (double)x->size[0]);
  emxFree_real_T(&x);
  return y;
}

/* End of code generation (rms.c) */
