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
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>

/* Function Definitions */
void gradient(const emxArray_real_T *x, emxArray_real_T *varargout_1)
{
  const double *x_data;
  double *varargout_1_data;
  int j;
  int loop_ub;
  int scalarLB;
  x_data = x->data;
  loop_ub = x->size[0];
  scalarLB = varargout_1->size[0];
  varargout_1->size[0] = x->size[0];
  emxEnsureCapacity_real_T(varargout_1, scalarLB);
  varargout_1_data = varargout_1->data;
  if (x->size[0] >= 2) {
    int vectorUB;
    varargout_1_data[0] = x_data[1] - x_data[0];
    scalarLB = (((x->size[0] - 2) / 2) << 1) + 2;
    vectorUB = scalarLB - 2;
    for (j = 2; j <= vectorUB; j += 2) {
      _mm_storeu_pd(&varargout_1_data[j - 1],
                    _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&x_data[j]),
                                          _mm_loadu_pd(&x_data[j - 2])),
                               _mm_set1_pd(2.0)));
    }
    for (j = scalarLB; j < loop_ub; j++) {
      varargout_1_data[j - 1] = (x_data[j] - x_data[j - 2]) / 2.0;
    }
    varargout_1_data[x->size[0] - 1] =
        x_data[x->size[0] - 1] - x_data[x->size[0] - 2];
  } else {
    scalarLB = varargout_1->size[0];
    varargout_1->size[0] = x->size[0];
    emxEnsureCapacity_real_T(varargout_1, scalarLB);
    varargout_1_data = varargout_1->data;
    for (j = 0; j < loop_ub; j++) {
      varargout_1_data[0] = 0.0;
    }
  }
}

/* End of code generation (gradient.c) */
