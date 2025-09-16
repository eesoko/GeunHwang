/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * periodogram.h
 *
 * Code generation for function 'periodogram'
 *
 */

#ifndef PERIODOGRAM_H
#define PERIODOGRAM_H

/* Include files */
#include "feature_extractor_codegen_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void periodogram(const double x_data[], int x_size, double varargin_3,
                 emxArray_real_T *Px, emxArray_real_T *w);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (periodogram.h) */
