/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CompactClassificationModel.c
 *
 * Code generation for function 'CompactClassificationModel'
 *
 */

/* Include files */
#include "CompactClassificationModel.h"
#include "cellstr_unique.h"
#include "convertCodes.h"
#include "predict_exercise_index_data.h"
#include "predict_exercise_index_internal_types.h"
#include "rt_nonfinite.h"
#include "strtrim.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void c_CompactClassificationModel_Co(c_classreg_learning_classif_Com *obj)
{
  static const char b_cv1[7] = {'P', 'u', 's', 'h', ' ', 'U', 'p'};
  static const char b_cv[5] = {'L', 'u', 'n', 'g', 'e'};
  static const char b_cv2[5] = {'S', 'q', 'u', 'a', 't'};
  cell_wrap_3 rv[6];
  cell_wrap_3 rv1[6];
  int i;
  unsigned char bcodes[6];
  rv[0].f1.size[0] = 1;
  rv[0].f1.size[1] = 13;
  for (i = 0; i < 13; i++) {
    rv[0].f1.data[i] = cv[i];
  }
  rv[1].f1.size[0] = 1;
  rv[1].f1.size[1] = 5;
  for (i = 0; i < 5; i++) {
    rv[1].f1.data[i] = b_cv[i];
  }
  rv[2].f1.size[0] = 1;
  rv[2].f1.size[1] = 14;
  for (i = 0; i < 14; i++) {
    rv[2].f1.data[i] = cv1[i];
  }
  rv[3].f1.size[0] = 1;
  rv[3].f1.size[1] = 7;
  for (i = 0; i < 7; i++) {
    rv[3].f1.data[i] = b_cv1[i];
  }
  rv[4].f1.size[0] = 1;
  rv[4].f1.size[1] = 18;
  for (i = 0; i < 18; i++) {
    rv[4].f1.data[i] = cv2[i];
  }
  rv[5].f1.size[0] = 1;
  rv[5].f1.size[1] = 5;
  for (i = 0; i < 5; i++) {
    rv[5].f1.data[i] = b_cv2[i];
  }
  rv1[0].f1.size[0] = 1;
  rv1[0].f1.size[1] = 13;
  for (i = 0; i < 13; i++) {
    rv1[0].f1.data[i] = cv[i];
  }
  rv1[1].f1.size[0] = 1;
  rv1[1].f1.size[1] = 5;
  for (i = 0; i < 5; i++) {
    rv1[1].f1.data[i] = b_cv[i];
  }
  rv1[2].f1.size[0] = 1;
  rv1[2].f1.size[1] = 14;
  for (i = 0; i < 14; i++) {
    rv1[2].f1.data[i] = cv1[i];
  }
  rv1[3].f1.size[0] = 1;
  rv1[3].f1.size[1] = 7;
  for (i = 0; i < 7; i++) {
    rv1[3].f1.data[i] = b_cv1[i];
  }
  rv1[4].f1.size[0] = 1;
  rv1[4].f1.size[1] = 18;
  for (i = 0; i < 18; i++) {
    rv1[4].f1.data[i] = cv2[i];
  }
  rv1[5].f1.size[0] = 1;
  rv1[5].f1.size[1] = 5;
  for (i = 0; i < 5; i++) {
    rv1[5].f1.data[i] = b_cv2[i];
  }
  categorical_convertCodes(rv1, rv, bcodes);
  for (i = 0; i < 6; i++) {
    if (bcodes[i] == 0) {
      bcodes[i] = MAX_uint8_T;
    }
  }
  for (i = 0; i < 6; i++) {
    int iRowS;
    boolean_T exitg1;
    obj->ClassLogicalIndices[i] = false;
    iRowS = 0;
    exitg1 = false;
    while ((!exitg1) && (iRowS < 6)) {
      if (i + 1 == bcodes[iRowS]) {
        obj->ClassLogicalIndices[i] = true;
        exitg1 = true;
      } else {
        iRowS++;
      }
    }
  }
}

int c_CompactClassificationModel_mi(const double obj_Cost[36],
                                    const double scoresIn_data[],
                                    const int scoresIn_size[2],
                                    double classnum_data[], double cost_data[],
                                    int cost_size[2])
{
  static const char b_cv1[7] = {'P', 'u', 's', 'h', ' ', 'U', 'p'};
  static const char b_cv[5] = {'L', 'u', 'n', 'g', 'e'};
  static const char b_cv2[5] = {'S', 'q', 'u', 'a', 't'};
  cell_wrap_3 names[6];
  cell_wrap_3 rv[6];
  emxArray_cell_wrap_3_6 r;
  double ib_data[6];
  double ex;
  int classnum_size;
  int cost_data_tmp;
  int i;
  int idx;
  int k;
  classnum_size = scoresIn_size[0];
  cost_size[0] = classnum_size;
  cost_size[1] = 6;
  for (k = 0; k < 6; k++) {
    if (classnum_size - 1 >= 0) {
      cost_data[classnum_size * k] = 0.0;
    }
    for (idx = 0; idx < 6; idx++) {
      ex = obj_Cost[idx + 6 * k];
      for (i = 0; i < classnum_size; i++) {
        cost_data_tmp = classnum_size * k;
        cost_data[cost_data_tmp] += scoresIn_data[scoresIn_size[0] * idx] * ex;
      }
    }
  }
  for (idx = 0; idx < classnum_size; idx++) {
    int b_idx;
    if (!rtIsNaN(cost_data[0])) {
      b_idx = 1;
    } else {
      boolean_T exitg1;
      b_idx = 0;
      cost_data_tmp = 2;
      exitg1 = false;
      while ((!exitg1) && (cost_data_tmp < 7)) {
        if (!rtIsNaN(cost_data[classnum_size * (cost_data_tmp - 1)])) {
          b_idx = cost_data_tmp;
          exitg1 = true;
        } else {
          cost_data_tmp++;
        }
      }
    }
    if (b_idx == 0) {
      b_idx = 1;
    } else {
      ex = cost_data[classnum_size * (b_idx - 1)];
      cost_data_tmp = b_idx + 1;
      for (k = cost_data_tmp; k < 7; k++) {
        double d;
        d = cost_data[classnum_size * (k - 1)];
        if (ex > d) {
          ex = d;
          b_idx = k;
        }
      }
    }
    classnum_data[0] = b_idx;
  }
  names[0].f1.size[0] = 1;
  names[0].f1.size[1] = 13;
  for (k = 0; k < 13; k++) {
    names[0].f1.data[k] = cv[k];
  }
  names[1].f1.size[0] = 1;
  names[1].f1.size[1] = 5;
  for (k = 0; k < 5; k++) {
    names[1].f1.data[k] = b_cv[k];
  }
  names[2].f1.size[0] = 1;
  names[2].f1.size[1] = 14;
  for (k = 0; k < 14; k++) {
    names[2].f1.data[k] = cv1[k];
  }
  names[3].f1.size[0] = 1;
  names[3].f1.size[1] = 7;
  for (k = 0; k < 7; k++) {
    names[3].f1.data[k] = b_cv1[k];
  }
  names[4].f1.size[0] = 1;
  names[4].f1.size[1] = 18;
  for (k = 0; k < 18; k++) {
    names[4].f1.data[k] = cv2[k];
  }
  names[5].f1.size[0] = 1;
  names[5].f1.size[1] = 5;
  for (k = 0; k < 5; k++) {
    names[5].f1.data[k] = b_cv2[k];
  }
  for (k = 0; k < 6; k++) {
    strtrim(names[k].f1.data, names[k].f1.size, rv[k].f1.data, rv[k].f1.size);
  }
  b_cellstr_unique(rv, r.data);
  cellstr_unique(rv, names, ib_data, &cost_data_tmp);
  return classnum_size;
}

/* End of code generation (CompactClassificationModel.c) */
