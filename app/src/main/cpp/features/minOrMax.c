/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * minOrMax.c
 *
 * Code generation for function 'minOrMax'
 *
 */

/* Include files */
#include "minOrMax.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"

/* Function Definitions */
double b_maximum(const emxArray_real_T *x, int *idx)
{
  const double *x_data;
  double ex;
  int b_k;
  int last;
  x_data = x->data;
  last = x->size[0];
  if (x->size[0] <= 2) {
    if (x->size[0] == 1) {
      ex = x_data[0];
      *idx = 1;
    } else {
      ex = x_data[x->size[0] - 1];
      if ((x_data[0] < ex) || (rtIsNaN(x_data[0]) && (!rtIsNaN(ex)))) {
        *idx = x->size[0];
      } else {
        ex = x_data[0];
        *idx = 1;
      }
    }
  } else {
    int b_idx;
    if (!rtIsNaN(x_data[0])) {
      b_idx = 1;
    } else {
      int k;
      boolean_T exitg1;
      b_idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!rtIsNaN(x_data[k - 1])) {
          b_idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (b_idx == 0) {
      ex = x_data[0];
      *idx = 1;
    } else {
      ex = x_data[b_idx - 1];
      *idx = b_idx;
      b_idx++;
      for (b_k = b_idx; b_k <= last; b_k++) {
        double d;
        d = x_data[b_k - 1];
        if (ex < d) {
          ex = d;
          *idx = b_k;
        }
      }
    }
  }
  return ex;
}

double maximum(const double x_data[], int x_size)
{
  double ex;
  int b_k;
  if (x_size <= 2) {
    if (x_size == 1) {
      ex = x_data[0];
    } else {
      ex = x_data[x_size - 1];
      if ((!(x_data[0] < ex)) && ((!rtIsNaN(x_data[0])) || rtIsNaN(ex))) {
        ex = x_data[0];
      }
    }
  } else {
    int idx;
    if (!rtIsNaN(x_data[0])) {
      idx = 1;
    } else {
      int k;
      boolean_T exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= x_size)) {
        if (!rtIsNaN(x_data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      ex = x_data[0];
    } else {
      ex = x_data[idx - 1];
      idx++;
      for (b_k = idx; b_k <= x_size; b_k++) {
        double d;
        d = x_data[b_k - 1];
        if (ex < d) {
          ex = d;
        }
      }
    }
  }
  return ex;
}

double minimum(const double x_data[], int x_size)
{
  double ex;
  int b_k;
  if (x_size <= 2) {
    if (x_size == 1) {
      ex = x_data[0];
    } else {
      ex = x_data[x_size - 1];
      if ((!(x_data[0] > ex)) && ((!rtIsNaN(x_data[0])) || rtIsNaN(ex))) {
        ex = x_data[0];
      }
    }
  } else {
    int idx;
    if (!rtIsNaN(x_data[0])) {
      idx = 1;
    } else {
      int k;
      boolean_T exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= x_size)) {
        if (!rtIsNaN(x_data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      ex = x_data[0];
    } else {
      ex = x_data[idx - 1];
      idx++;
      for (b_k = idx; b_k <= x_size; b_k++) {
        double d;
        d = x_data[b_k - 1];
        if (ex > d) {
          ex = d;
        }
      }
    }
  }
  return ex;
}

/* End of code generation (minOrMax.c) */
