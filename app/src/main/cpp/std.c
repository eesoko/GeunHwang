/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * std.c
 *
 * Code generation for function 'std'
 *
 */

/* Include files */
#include "std.h"
#include "blockedSummation.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
double b_std(const double x_data[], int x_size)
{
  emxArray_real_T b_x_data;
  double y;
  int k;
  if (x_size == 0) {
    y = rtNaN;
  } else if (x_size == 1) {
    if ((!rtIsInf(x_data[0])) && (!rtIsNaN(x_data[0]))) {
      y = 0.0;
    } else {
      y = rtNaN;
    }
  } else {
    double scale;
    double xbar;
    b_x_data.data = (double *)&x_data[0];
    b_x_data.size = &x_size;
    b_x_data.allocatedSize = -1;
    b_x_data.numDimensions = 1;
    b_x_data.canFreeData = false;
    xbar = colMajorFlatIter(&b_x_data, x_size) / (double)x_size;
    y = 0.0;
    scale = 3.3121686421112381E-170;
    for (k = 0; k < x_size; k++) {
      double d;
      d = fabs(x_data[k] - xbar);
      if (d > scale) {
        double t;
        t = scale / d;
        y = y * t * t + 1.0;
        scale = d;
      } else {
        double t;
        t = d / scale;
        y += t * t;
      }
    }
    y = scale * sqrt(y);
    y /= sqrt((double)x_size - 1.0);
  }
  return y;
}

/* End of code generation (std.c) */
