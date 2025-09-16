/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * predict_exercise.h
 *
 * Code generation for function 'predict_exercise'
 *
 */

#ifndef PREDICT_EXERCISE_H
#define PREDICT_EXERCISE_H

/* Include files */
#include "predict_exercise_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void predict_exercise(const double input_features[32],
                             categorical *label);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (predict_exercise.h) */
