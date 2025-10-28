/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * pca.h
 *
 * Code generation for function 'pca'
 *
 */

#ifndef PCA_H
#define PCA_H

/* Include files */
#include "feature_extractor_codegen_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
int pca(const emxArray_real_T *x, double varargout_1_data[],
        int varargout_1_size[2], emxArray_real_T *varargout_2,
        double varargout_3_data[]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (pca.h) */
