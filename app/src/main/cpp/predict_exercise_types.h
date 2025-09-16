/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * predict_exercise_types.h
 *
 * Code generation for function 'predict_exercise'
 *
 */

#ifndef PREDICT_EXERCISE_TYPES_H
#define PREDICT_EXERCISE_TYPES_H

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef struct_emxArray_char_T_1x18
#define struct_emxArray_char_T_1x18
struct emxArray_char_T_1x18 {
  char data[18];
  int size[2];
};
#endif /* struct_emxArray_char_T_1x18 */
#ifndef typedef_emxArray_char_T_1x18
#define typedef_emxArray_char_T_1x18
typedef struct emxArray_char_T_1x18 emxArray_char_T_1x18;
#endif /* typedef_emxArray_char_T_1x18 */

#ifndef typedef_cell_wrap_0
#define typedef_cell_wrap_0
typedef struct {
  emxArray_char_T_1x18 f1;
} cell_wrap_0;
#endif /* typedef_cell_wrap_0 */

#ifndef typedef_categorical
#define typedef_categorical
typedef struct {
  unsigned char codes;
  cell_wrap_0 categoryNames[6];
} categorical;
#endif /* typedef_categorical */

#endif
/* End of code generation (predict_exercise_types.h) */
