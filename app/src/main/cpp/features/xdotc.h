/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xdotc.h
 *
 * Code generation for function 'xdotc'
 *
 */

#ifndef XDOTC_H
#define XDOTC_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
double b_xdotc(const double x[9], const double y[9], int iy0);

double xdotc(int n, const double x_data[], int ix0, const double y_data[],
             int iy0);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (xdotc.h) */
