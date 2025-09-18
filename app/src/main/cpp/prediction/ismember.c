/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ismember.c
 *
 * Code generation for function 'ismember'
 *
 */

/* Include files */
#include "ismember.h"
#include "convertCodes.h"
#include "predict_exercise_data.h"
#include "predict_exercise_internal_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void categorical_ismember(boolean_T tf[6], double loc[6])
{
  static const char b_cv1[7] = {'P', 'u', 's', 'h', ' ', 'U', 'p'};
  static const char b_cv[5] = {'L', 'u', 'n', 'g', 'e'};
  static const char b_cv2[5] = {'S', 'q', 'u', 'a', 't'};
  cell_wrap_3 rv[6];
  cell_wrap_3 rv1[6];
  int b_loc[6];
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
    tf[i] = false;
    b_loc[i] = 0;
    iRowS = 0;
    exitg1 = false;
    while ((!exitg1) && (iRowS < 6)) {
      if (i + 1 == bcodes[iRowS]) {
        tf[i] = true;
        b_loc[i] = iRowS + 1;
        exitg1 = true;
      } else {
        iRowS++;
      }
    }
    loc[i] = b_loc[i];
  }
}

/* End of code generation (ismember.c) */
