/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sqrt.c
 *
 * Code generation for function 'sqrt'
 *
 */

/* Include files */
#include "sqrt.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>

/* Function Definitions */
void b_sqrt(double x_data[], const int *x_size)
{
  int i;
  int k;
  int scalarLB;
  int vectorUB;
  i = *x_size;
  scalarLB = (i / 2) << 1;
  vectorUB = scalarLB - 2;
  for (k = 0; k <= vectorUB; k += 2) {
    __m128d r;
    r = _mm_loadu_pd(&x_data[k]);
    _mm_storeu_pd(&x_data[k], _mm_sqrt_pd(r));
  }
  for (k = scalarLB; k < i; k++) {
    x_data[k] = sqrt(x_data[k]);
  }
}

/* End of code generation (sqrt.c) */
