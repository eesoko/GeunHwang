/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * feature_extractor_codegen_initialize.c
 *
 * Code generation for function 'feature_extractor_codegen_initialize'
 *
 */

/* Include files */
#include "feature_extractor_codegen_initialize.h"
#include "feature_extractor_codegen_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
void feature_extractor_codegen_initialize(void)
{
  omp_init_nest_lock(&feature_extractor_codegen_nestLockGlobal);
  isInitialized_feature_extractor_codegen = true;
}

/* End of code generation (feature_extractor_codegen_initialize.c) */
