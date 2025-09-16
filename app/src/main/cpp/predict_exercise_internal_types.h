/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * predict_exercise_internal_types.h
 *
 * Code generation for function 'predict_exercise'
 *
 */

#ifndef PREDICT_EXERCISE_INTERNAL_TYPES_H
#define PREDICT_EXERCISE_INTERNAL_TYPES_H

/* Include files */
#include "predict_exercise_types.h"
#include "rtwtypes.h"

/* Type Definitions */
#ifndef c_enum_c_classreg_learning_code
#define c_enum_c_classreg_learning_code
enum c_classreg_learning_coderutils_
{
  linear = 1, /* Default value */
  gaussian = 2,
  rbf = 2,
  polynomial = 3
};
#endif /* c_enum_c_classreg_learning_code */
#ifndef c_typedef_c_classreg_learning_c
#define c_typedef_c_classreg_learning_c
typedef enum c_classreg_learning_coderutils_ c_classreg_learning_coderutils_;
#endif /* c_typedef_c_classreg_learning_c */

#endif
/* End of code generation (predict_exercise_internal_types.h) */
