/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * convertCodes.h
 *
 * Code generation for function 'convertCodes'
 *
 */

#ifndef CONVERTCODES_H
#define CONVERTCODES_H

/* Include files */
#include "predict_exercise_internal_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void categorical_convertCodes(const cell_wrap_3 bnames[6],
                              const cell_wrap_3 anames[6],
                              unsigned char bcodesout[6]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (convertCodes.h) */
