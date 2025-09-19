/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * predict_exercise_index_internal_types.h
 *
 * Code generation for function 'predict_exercise_index'
 *
 */

#ifndef PREDICT_EXERCISE_INDEX_INTERNAL_TYPES_H
#define PREDICT_EXERCISE_INDEX_INTERNAL_TYPES_H

/* Include files */
#include "predict_exercise_index_types.h"
#include "rtwtypes.h"

/* Type Definitions */
#ifndef c_enum_c_classreg_learning_code
#define c_enum_c_classreg_learning_code
enum c_classreg_learning_coderutils_
{
  Logit = 0, /* Default value */
  Doublelogit,
  Invlogit,
  Ismax,
  Sign,
  Symmetric,
  Symmetricismax,
  Symmetriclogit,
  Identity
};
#endif /* c_enum_c_classreg_learning_code */
#ifndef c_typedef_c_classreg_learning_c
#define c_typedef_c_classreg_learning_c
typedef enum c_classreg_learning_coderutils_ c_classreg_learning_coderutils_;
#endif /* c_typedef_c_classreg_learning_c */

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

#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3
typedef struct {
  emxArray_char_T_1x18 f1;
} cell_wrap_3;
#endif /* typedef_cell_wrap_3 */

#ifndef typedef_emxArray_cell_wrap_3_6
#define typedef_emxArray_cell_wrap_3_6
typedef struct {
  cell_wrap_3 data[6];
} emxArray_cell_wrap_3_6;
#endif /* typedef_emxArray_cell_wrap_3_6 */

#ifndef d_typedef_c_classreg_learning_c
#define d_typedef_c_classreg_learning_c
typedef struct {
  double Prior[6];
  boolean_T ClassLogicalIndices[6];
  double Cost[36];
  double LogDetSigma;
} c_classreg_learning_classif_Com;
#endif /* d_typedef_c_classreg_learning_c */

#endif
/* End of code generation (predict_exercise_index_internal_types.h) */
