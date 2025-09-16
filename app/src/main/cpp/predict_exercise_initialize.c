/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * predict_exercise_initialize.c
 *
 * Code generation for function 'predict_exercise_initialize'
 *
 */

/* Include files */
#include "predict_exercise_initialize.h"
#include "predict_exercise_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
void predict_exercise_initialize(void)
{
  omp_init_nest_lock(&predict_exercise_nestLockGlobal);
  isInitialized_predict_exercise = true;
}

/* End of code generation (predict_exercise_initialize.c) */
