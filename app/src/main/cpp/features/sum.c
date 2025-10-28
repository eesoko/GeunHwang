/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "sum.h"
#include "blockedSummation.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>

/* Function Definitions */
double b_sum(const double x_data[], int x_size)
{
  emxArray_real_T b_x_data;
  double y;
  if (x_size == 0) {
    y = 0.0;
  } else {
    b_x_data.data = (double *)&x_data[0];
    b_x_data.size = &x_size;
    b_x_data.allocatedSize = -1;
    b_x_data.numDimensions = 1;
    b_x_data.canFreeData = false;
    y = colMajorFlatIter(&b_x_data, x_size);
  }
  return y;
}

void sum(const emxArray_real_T *x, emxArray_real_T *y)
{
  const double *x_data;
  double *y_data;
  int b_xj;
  int xj;
  x_data = x->data;
  if (x->size[0] == 0) {
    y->size[0] = 0;
  } else {
    int scalarLB;
    int vectorUB;
    int vstride;
    vstride = x->size[0];
    scalarLB = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real_T(y, scalarLB);
    y_data = y->data;
    for (xj = 0; xj < vstride; xj++) {
      y_data[xj] = x_data[xj];
    }
    scalarLB = (x->size[0] / 2) << 1;
    vectorUB = scalarLB - 2;
    for (xj = 0; xj < 2; xj++) {
      int xoffset;
      xoffset = (xj + 1) * vstride;
      for (b_xj = 0; b_xj <= vectorUB; b_xj += 2) {
        __m128d r;
        r = _mm_loadu_pd(&y_data[b_xj]);
        _mm_storeu_pd(&y_data[b_xj],
                      _mm_add_pd(r, _mm_loadu_pd(&x_data[xoffset + b_xj])));
      }
      for (b_xj = scalarLB; b_xj < vstride; b_xj++) {
        y_data[b_xj] += x_data[xoffset + b_xj];
      }
    }
  }
}

/* End of code generation (sum.c) */
