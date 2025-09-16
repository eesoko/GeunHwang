/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * loss.c
 *
 * Code generation for function 'loss'
 *
 */

/* Include files */
#include "loss.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void loss(const double M[90], const double pscore[15], double vloss[6])
{
  int b_k;
  int k;
  int p;
  boolean_T exitg1;
  boolean_T valueisnan;
  valueisnan = true;
  p = 0;
  exitg1 = false;
  while ((!exitg1) && (p < 15)) {
    if (!rtIsNaN(pscore[p])) {
      valueisnan = false;
      exitg1 = true;
    } else {
      p++;
    }
  }
  for (k = 0; k < 6; k++) {
    vloss[k] = rtNaN;
  }
  if (!valueisnan) {
    for (b_k = 0; b_k < 6; b_k++) {
      double maxval[15];
      double y;
      for (k = 0; k < 15; k++) {
        y = 1.0 - M[b_k + 6 * k] * pscore[k];
        if (rtIsNaN(y) || (y > 0.0)) {
          maxval[k] = y;
        } else {
          maxval[k] = 0.0;
        }
      }
      if (rtIsNaN(maxval[0])) {
        y = 0.0;
        p = 0;
      } else {
        y = maxval[0];
        p = 1;
      }
      for (k = 0; k < 14; k++) {
        double d;
        d = maxval[k + 1];
        if (!rtIsNaN(d)) {
          y += d;
          p++;
        }
      }
      vloss[b_k] = y / (double)p / 2.0;
    }
  }
}

/* End of code generation (loss.c) */
