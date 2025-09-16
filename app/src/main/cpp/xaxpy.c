/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xaxpy.c
 *
 * Code generation for function 'xaxpy'
 *
 */

/* Include files */
#include "xaxpy.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>

/* Function Definitions */
void b_xaxpy(double a, double y[9], int iy0)
{
  if (!(a == 0.0)) {
    y[iy0 - 1] += a * y[1];
    y[iy0] += a * y[2];
  }
}

void c_xaxpy(int n, double a, const double x_data[], int ix0, double y_data[])
{
  int k;
  if ((n >= 1) && (!(a == 0.0))) {
    int scalarLB;
    int vectorUB;
    scalarLB = (n / 2) << 1;
    vectorUB = scalarLB - 2;
    for (k = 0; k <= vectorUB; k += 2) {
      __m128d r;
      r = _mm_loadu_pd(&y_data[k + 1]);
      _mm_storeu_pd(
          &y_data[k + 1],
          _mm_add_pd(r, _mm_mul_pd(_mm_set1_pd(a),
                                   _mm_loadu_pd(&x_data[(ix0 + k) - 1]))));
    }
    for (k = scalarLB; k < n; k++) {
      y_data[k + 1] += a * x_data[(ix0 + k) - 1];
    }
  }
}

void d_xaxpy(int n, double a, const double x_data[], double y_data[], int iy0)
{
  int k;
  if ((n >= 1) && (!(a == 0.0))) {
    int scalarLB;
    int vectorUB;
    scalarLB = (n / 2) << 1;
    vectorUB = scalarLB - 2;
    for (k = 0; k <= vectorUB; k += 2) {
      __m128d r;
      int i;
      i = (iy0 + k) - 1;
      r = _mm_loadu_pd(&y_data[i]);
      _mm_storeu_pd(&y_data[i],
                    _mm_add_pd(r, _mm_mul_pd(_mm_set1_pd(a),
                                             _mm_loadu_pd(&x_data[k + 1]))));
    }
    for (k = scalarLB; k < n; k++) {
      vectorUB = (iy0 + k) - 1;
      y_data[vectorUB] += a * x_data[k + 1];
    }
  }
}

void xaxpy(int n, double a, int ix0, double y_data[], int iy0)
{
  int k;
  if ((n >= 1) && (!(a == 0.0))) {
    for (k = 0; k < n; k++) {
      int i;
      i = (iy0 + k) - 1;
      y_data[i] += a * y_data[(ix0 + k) - 1];
    }
  }
}

/* End of code generation (xaxpy.c) */
