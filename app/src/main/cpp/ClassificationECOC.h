/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ClassificationECOC.h
 *
 * Code generation for function 'ClassificationECOC'
 *
 */

#ifndef CLASSIFICATIONECOC_H
#define CLASSIFICATIONECOC_H

/* Include files */
#include "predict_exercise_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
unsigned char ClassificationECOC_predict(const double Xin[32],
                                         cell_wrap_0 labels_categoryNames[6]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (ClassificationECOC.h) */
