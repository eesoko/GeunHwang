/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gradient.c
 *
 * Code generation for function 'gradient'
 *
 */

/* Include files */
#include "gradient.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>

/* Function Definitions */
int gradient(const double x_data[], int x_size, double varargout_1_data[])
{
  int j;
  int varargout_1_size;
  varargout_1_size = x_size;
  if (x_size >= 2) {
    int scalarLB;
    int vectorUB;
    varargout_1_data[0] = x_data[1] - x_data[0];
    scalarLB = (((x_size - 2) / 2) << 1) + 2;
    vectorUB = scalarLB - 2;
    for (j = 2; j <= vectorUB; j += 2) {
      _mm_storeu_pd(&varargout_1_data[j - 1],
                    _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&x_data[j]),
                                          _mm_loadu_pd(&x_data[j - 2])),
                               _mm_set1_pd(2.0)));
    }
    for (j = scalarLB; j < x_size; j++) {
      varargout_1_data[j - 1] = (x_data[j] - x_data[j - 2]) / 2.0;
    }
    varargout_1_data[x_size - 1] = x_data[x_size - 1] - x_data[x_size - 2];
  } else if (x_size - 1 >= 0) {
    varargout_1_data[0] = 0.0;
  }
  return varargout_1_size;
}

/* End of code generation (gradient.c) */
