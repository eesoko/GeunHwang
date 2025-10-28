/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlascl.h
 *
 * Code generation for function 'xzlascl'
 *
 */

#ifndef XZLASCL_H
#define XZLASCL_H

/* Include files */
#include "feature_extractor_codegen_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_xzlascl(double cfrom, double cto, int m, emxArray_real_T *A, int lda);

void xzlascl(double cfrom, double cto, int m, double A_data[]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (xzlascl.h) */
