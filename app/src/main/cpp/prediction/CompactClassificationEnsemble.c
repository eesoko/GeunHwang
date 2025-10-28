/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CompactClassificationEnsemble.c
 *
 * Code generation for function 'CompactClassificationEnsemble'
 *
 */

/* Include files */
#include "CompactClassificationEnsemble.h"
#include "CompactEnsemble.h"
#include "cellstr_unique.h"
#include "predict_exercise_index_data.h"
#include "predict_exercise_index_internal_types.h"
#include "rt_nonfinite.h"
#include "strtrim.h"

/* Function Definitions */
void c_CompactClassificationEnsemble(const double Xin[32], double score[6])
{
  static const char b_cv1[7] = {'P', 'u', 's', 'h', ' ', 'U', 'p'};
  static const char b_cv[5] = {'L', 'u', 'n', 'g', 'e'};
  static const char b_cv2[5] = {'S', 'q', 'u', 'a', 't'};
  cell_wrap_3 names[6];
  cell_wrap_3 rv[6];
  emxArray_cell_wrap_3_6 r;
  int i;
  CompactEnsemble_ensemblePredict(Xin, score);
  names[0].f1.size[0] = 1;
  names[0].f1.size[1] = 13;
  for (i = 0; i < 13; i++) {
    names[0].f1.data[i] = cv[i];
  }
  names[1].f1.size[0] = 1;
  names[1].f1.size[1] = 5;
  for (i = 0; i < 5; i++) {
    names[1].f1.data[i] = b_cv[i];
  }
  names[2].f1.size[0] = 1;
  names[2].f1.size[1] = 14;
  for (i = 0; i < 14; i++) {
    names[2].f1.data[i] = cv1[i];
  }
  names[3].f1.size[0] = 1;
  names[3].f1.size[1] = 7;
  for (i = 0; i < 7; i++) {
    names[3].f1.data[i] = b_cv1[i];
  }
  names[4].f1.size[0] = 1;
  names[4].f1.size[1] = 18;
  for (i = 0; i < 18; i++) {
    names[4].f1.data[i] = cv2[i];
  }
  names[5].f1.size[0] = 1;
  names[5].f1.size[1] = 5;
  for (i = 0; i < 5; i++) {
    names[5].f1.data[i] = b_cv2[i];
  }
  for (i = 0; i < 6; i++) {
    strtrim(names[i].f1.data, names[i].f1.size, rv[i].f1.data, rv[i].f1.size);
  }
  b_cellstr_unique(rv, r.data);
}

/* End of code generation (CompactClassificationEnsemble.c) */
