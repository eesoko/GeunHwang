/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlangeM.c
 *
 * Code generation for function 'xzlangeM'
 *
 */

/* Include files */
#include "xzlangeM.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
double xzlangeM(const emxArray_real_T *x)
{
  const double *x_data;
  double y;
  boolean_T b;
  x_data = x->data;
  y = 0.0;
  b = (x->size[0] == 0);
  if (!b) {
    int k;
    k = 0;
    int exitg1;
    int i;
    do {
      exitg1 = 0;
      i = x->size[0] * 3;
      if (k <= i - 1) {
        double absxk;
        absxk = fabs(x_data[k]);
        if (rtIsNaN(absxk)) {
          y = rtNaN;
          exitg1 = 1;
        } else {
          if (absxk > y) {
            y = absxk;
          }
          k++;
        }
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  return y;
}

/* End of code generation (xzlangeM.c) */
