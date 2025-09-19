/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * feature_extractor_codegen.h
 *
 * Code generation for function 'feature_extractor_codegen'
 *
 */

#ifndef FEATURE_EXTRACTOR_CODEGEN_H
#define FEATURE_EXTRACTOR_CODEGEN_H

/* Include files */
#include "feature_extractor_codegen_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void feature_extractor_codegen(const emxArray_real_T *raw_data,
                                      double Fs_actual, double features[32]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (feature_extractor_codegen.h) */
