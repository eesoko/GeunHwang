/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xaxpy.h
 *
 * Code generation for function 'xaxpy'
 *
 */

#ifndef XAXPY_H
#define XAXPY_H

/* Include files */
#include "feature_extractor_codegen_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_xaxpy(double a, double y[9], int iy0);

void c_xaxpy(int n, double a, const emxArray_real_T *x, int ix0,
             emxArray_real_T *y);

void d_xaxpy(int n, double a, const emxArray_real_T *x, emxArray_real_T *y,
             int iy0);

void xaxpy(int n, double a, int ix0, emxArray_real_T *y, int iy0);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (xaxpy.h) */
