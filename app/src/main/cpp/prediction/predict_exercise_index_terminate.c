/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * predict_exercise_index_terminate.c
 *
 * Code generation for function 'predict_exercise_index_terminate'
 *
 */

/* Include files */
#include "predict_exercise_index_terminate.h"
#include "predict_exercise_index_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
void predict_exercise_index_terminate(void)
{
  omp_destroy_nest_lock(&predict_exercise_index_nestLockGlobal);
  isInitialized_predict_exercise_index = false;
}

/* End of code generation (predict_exercise_index_terminate.c) */
