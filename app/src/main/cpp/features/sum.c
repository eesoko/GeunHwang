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
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <string.h>

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

int sum(const double x_data[], const int x_size[2], double y_data[])
{
  int k;
  int xj;
  int y_size;
  if (x_size[0] == 0) {
    y_size = 0;
  } else {
    int scalarLB;
    int vectorUB;
    y_size = x_size[0];
    memcpy(&y_data[0], &x_data[0], (unsigned int)y_size * sizeof(double));
    scalarLB = (x_size[0] / 2) << 1;
    vectorUB = scalarLB - 2;
    for (k = 0; k < 2; k++) {
      int xoffset;
      xoffset = (k + 1) * y_size;
      for (xj = 0; xj <= vectorUB; xj += 2) {
        __m128d r;
        r = _mm_loadu_pd(&y_data[xj]);
        _mm_storeu_pd(&y_data[xj],
                      _mm_add_pd(r, _mm_loadu_pd(&x_data[xoffset + xj])));
      }
      for (xj = scalarLB; xj < y_size; xj++) {
        y_data[xj] += x_data[xoffset + xj];
      }
    }
  }
  return y_size;
}

/* End of code generation (sum.c) */
