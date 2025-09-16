/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Poly.c
 *
 * Code generation for function 'Poly'
 *
 */

/* Include files */
#include "Poly.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>

/* Function Definitions */
void Poly(const double svT_data[], const int svT_size[2], const double x[32],
          double kernelProduct_data[], int kernelProduct_size[2])
{
  double d;
  int i;
  int i1;
  int loop_ub;
  int scalarLB;
  int vectorUB;
  loop_ub = svT_size[1];
  kernelProduct_size[1] = svT_size[1];
  for (i = 0; i < loop_ub; i++) {
    d = 0.0;
    for (i1 = 0; i1 < 32; i1++) {
      d += x[i1] * svT_data[i1 + 32 * i];
    }
    kernelProduct_data[i] = d + 1.0;
  }
  kernelProduct_size[0] = 1;
  loop_ub = svT_size[1] - 1;
  scalarLB = (svT_size[1] / 2) << 1;
  vectorUB = scalarLB - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    __m128d r;
    r = _mm_loadu_pd(&kernelProduct_data[i]);
    _mm_storeu_pd(&kernelProduct_data[i], _mm_mul_pd(r, r));
  }
  for (i = scalarLB; i <= loop_ub; i++) {
    d = kernelProduct_data[i];
    d *= d;
    kernelProduct_data[i] = d;
  }
}

/* End of code generation (Poly.c) */
