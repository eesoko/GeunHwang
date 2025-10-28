/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xdotc.c
 *
 * Code generation for function 'xdotc'
 *
 */

/* Include files */
#include "xdotc.h"
#include "rt_nonfinite.h"

/* Function Definitions */
double xdotc(const double x[9], const double y[9], int iy0)
{
  return x[1] * y[iy0 - 1] + x[2] * y[iy0];
}

/* End of code generation (xdotc.c) */
