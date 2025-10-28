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
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>

/* Function Definitions */
void corrcoef(const emxArray_real_T *x, const emxArray_real_T *varargin_1,
              double r[4])
{
  emxArray_real_T *result;
  double d[2];
  const double *varargin_1_data;
  const double *x_data;
  double absrij;
  double s;
  double *result_data;
  int b_j;
  int b_loop_ub;
  int i;
  int j;
  int loop_ub;
  int vectorUB;
  varargin_1_data = varargin_1->data;
  x_data = x->data;
  emxInit_real_T(&result, 2);
  loop_ub = x->size[0];
  b_loop_ub = result->size[0] * result->size[1];
  result->size[0] = x->size[0];
  result->size[1] = 2;
  emxEnsureCapacity_real_T(result, b_loop_ub);
  result_data = result->data;
  for (i = 0; i < loop_ub; i++) {
    result_data[i] = x_data[i];
  }
  b_loop_ub = varargin_1->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    result_data[i + result->size[0]] = varargin_1_data[i];
  }
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  r[3] = 0.0;
  if (result->size[0] < 2) {
    r[0] = rtNaN;
    r[1] = rtNaN;
    r[2] = rtNaN;
    r[3] = rtNaN;
  } else {
    int fm;
    b_loop_ub = (result->size[0] / 2) << 1;
    vectorUB = b_loop_ub - 2;
    for (j = 0; j < 2; j++) {
      s = 0.0;
      for (i = 0; i < loop_ub; i++) {
        s += result_data[i + result->size[0] * j];
      }
      s /= (double)result->size[0];
      for (i = 0; i <= vectorUB; i += 2) {
        __m128d b_r;
        b_r = _mm_loadu_pd(&result_data[i + result->size[0] * j]);
        _mm_storeu_pd(&result_data[i + result->size[0] * j],
                      _mm_sub_pd(b_r, _mm_set1_pd(s)));
      }
      for (i = b_loop_ub; i < loop_ub; i++) {
        result_data[i + result->size[0] * j] -= s;
      }
    }
    fm = result->size[0] - 1;
    for (b_j = 0; b_j < 2; b_j++) {
      s = 0.0;
      for (i = 0; i < loop_ub; i++) {
        absrij = result_data[i + result->size[0] * b_j];
        s += absrij * absrij;
      }
      b_loop_ub = b_j << 1;
      r[b_j + b_loop_ub] = s / (double)fm;
      vectorUB = b_j + 2;
      for (i = vectorUB; i < 3; i++) {
        s = 0.0;
        for (j = 0; j < loop_ub; j++) {
          s += result_data[j + result->size[0]] *
               result_data[j + result->size[0] * b_j];
        }
        r[b_loop_ub + 1] = s / (double)fm;
      }
    }
  }
  emxFree_real_T(&result);
  d[0] = sqrt(r[0]);
  d[1] = sqrt(r[3]);
  for (j = 0; j < 2; j++) {
    vectorUB = j + 2;
    for (i = vectorUB; i < 3; i++) {
      b_loop_ub = (j << 1) + 1;
      r[b_loop_ub] = r[b_loop_ub] / d[1] / d[j];
    }
    for (i = vectorUB; i < 3; i++) {
      b_loop_ub = (j << 1) + 1;
      s = r[b_loop_ub];
      absrij = fabs(s);
      if (absrij > 1.0) {
        s /= absrij;
        r[b_loop_ub] = s;
      }
      r[j + 2] = r[b_loop_ub];
    }
    b_loop_ub = j + (j << 1);
    s = r[b_loop_ub];
    if (s > 0.0) {
      if (rtIsNaN(s)) {
        r[b_loop_ub] = rtNaN;
      } else if (s < 0.0) {
        r[b_loop_ub] = -1.0;
      } else {
        r[b_loop_ub] = (s > 0.0);
      }
    } else {
      r[b_loop_ub] = rtNaN;
    }
  }
}

/* End of code generation (corrcoef.c) */
