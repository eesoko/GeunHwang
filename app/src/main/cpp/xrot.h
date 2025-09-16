/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xrot.h
 *
 * Code generation for function 'xrot'
 *
 */

#ifndef XROT_H
#define XROT_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_xrot(int n, double x_data[], int ix0, int iy0, double c, double s);

void xrot(double x[9], int ix0, int iy0, double c, double s);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (xrot.h) */
