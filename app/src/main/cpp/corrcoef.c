/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * corrcoef.c
 *
 * Code generation for function 'corrcoef'
 *
 */

/* Include files */
#include "corrcoef.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Definitions */
void corrcoef(const double x_data[], int x_size, const double varargin_1_data[],
              int varargin_1_size, double r[4])
{
  double result_data[300];
  double d[2];
  double absrij;
  double s;
  int b_j;
  int i;
  int j;
  int r_tmp;
  int scalarLB;
  if (x_size - 1 >= 0) {
    memcpy(&result_data[0], &x_data[0], (unsigned int)x_size * sizeof(double));
  }
  for (i = 0; i < varargin_1_size; i++) {
    result_data[i + x_size] = varargin_1_data[i];
  }
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  r[3] = 0.0;
  if (x_size < 2) {
    r[0] = rtNaN;
    r[1] = rtNaN;
    r[2] = rtNaN;
    r[3] = rtNaN;
  } else {
    int vectorUB;
    scalarLB = (x_size / 2) << 1;
    vectorUB = scalarLB - 2;
    for (j = 0; j < 2; j++) {
      s = 0.0;
      for (i = 0; i < x_size; i++) {
        s += result_data[i + x_size * j];
      }
      s /= (double)x_size;
      for (i = 0; i <= vectorUB; i += 2) {
        __m128d b_r;
        r_tmp = i + x_size * j;
        b_r = _mm_loadu_pd(&result_data[r_tmp]);
        _mm_storeu_pd(&result_data[r_tmp], _mm_sub_pd(b_r, _mm_set1_pd(s)));
      }
      for (i = scalarLB; i < x_size; i++) {
        r_tmp = i + x_size * j;
        result_data[r_tmp] -= s;
      }
    }
    for (j = 0; j < 2; j++) {
      s = 0.0;
      for (i = 0; i < x_size; i++) {
        absrij = result_data[i + x_size * j];
        s += absrij * absrij;
      }
      r_tmp = j << 1;
      r[j + r_tmp] = s / (double)(x_size - 1);
      scalarLB = j + 2;
      for (i = scalarLB; i < 3; i++) {
        s = 0.0;
        for (b_j = 0; b_j < x_size; b_j++) {
          s += result_data[b_j + x_size] * result_data[b_j + x_size * j];
        }
        r[r_tmp + 1] = s / (double)(x_size - 1);
      }
    }
  }
  d[0] = sqrt(r[0]);
  d[1] = sqrt(r[3]);
  for (b_j = 0; b_j < 2; b_j++) {
    scalarLB = b_j + 2;
    for (i = scalarLB; i < 3; i++) {
      r_tmp = (b_j << 1) + 1;
      r[r_tmp] = r[r_tmp] / d[1] / d[b_j];
    }
    for (i = scalarLB; i < 3; i++) {
      r_tmp = (b_j << 1) + 1;
      s = r[r_tmp];
      absrij = fabs(s);
      if (absrij > 1.0) {
        s /= absrij;
        r[r_tmp] = s;
      }
      r[b_j + 2] = r[r_tmp];
    }
    r_tmp = b_j + (b_j << 1);
    s = r[r_tmp];
    if (s > 0.0) {
      if (rtIsNaN(s)) {
        r[r_tmp] = rtNaN;
      } else if (s < 0.0) {
        r[r_tmp] = -1.0;
      } else {
        r[r_tmp] = (s > 0.0);
      }
    } else {
      r[r_tmp] = rtNaN;
    }
  }
}

/* End of code generation (corrcoef.c) */
