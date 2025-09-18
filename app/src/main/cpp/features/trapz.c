/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * trapz.c
 *
 * Code generation for function 'trapz'
 *
 */

/* Include files */
#include "trapz.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>

/* Function Definitions */
double trapz(const double x_data[], int x_size, const double y_data[],
             int y_size)
{
  double c_data[150];
  double z;
  int ia;
  int iac;
  z = 0.0;
  if (y_size <= 1) {
    if ((x_size == 1) && (rtIsInf(x_data[0]) || rtIsNaN(x_data[0]))) {
      z = rtNaN;
    }
  } else {
    double c1;
    if (x_size == 1) {
      for (iac = 0; iac < y_size; iac++) {
        c_data[iac] = x_data[0];
      }
      c1 = 0.5 * x_data[0];
      c_data[0] = c1;
      c_data[y_size - 1] = c1;
    } else {
      int scalarLB;
      int vectorUB;
      c_data[0] = 0.5 * (x_data[1] - x_data[0]);
      scalarLB = (((y_size - 2) / 2) << 1) + 2;
      vectorUB = scalarLB - 2;
      for (iac = 2; iac <= vectorUB; iac += 2) {
        _mm_storeu_pd(&c_data[iac - 1],
                      _mm_mul_pd(_mm_set1_pd(0.5),
                                 _mm_sub_pd(_mm_loadu_pd(&x_data[iac]),
                                            _mm_loadu_pd(&x_data[iac - 2]))));
      }
      for (iac = scalarLB; iac < y_size; iac++) {
        c_data[iac - 1] = 0.5 * (x_data[iac] - x_data[iac - 2]);
      }
      c_data[y_size - 1] = 0.5 * (x_data[y_size - 1] - x_data[y_size - 2]);
    }
    for (iac = 1; y_size < 0 ? iac >= 1 : iac <= 1; iac += y_size) {
      c1 = 0.0;
      for (ia = 1; ia <= y_size; ia++) {
        c1 += y_data[ia - 1] * c_data[ia - 1];
      }
      z += c1;
    }
  }
  return z;
}

/* End of code generation (trapz.c) */
