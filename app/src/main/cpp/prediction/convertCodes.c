/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * convertCodes.c
 *
 * Code generation for function 'convertCodes'
 *
 */

/* Include files */
#include "convertCodes.h"
#include "cellstr_unique.h"
#include "find.h"
#include "introsort.h"
#include "predict_exercise_internal_types.h"
#include "rt_nonfinite.h"
#include "strcmp.h"
#include "string1.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Declarations */
static double rt_roundd_snf(double u);

/* Function Definitions */
static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }
  return y;
}

void categorical_convertCodes(const cell_wrap_3 bnames[6],
                              const cell_wrap_3 anames[6],
                              unsigned char bcodesout[6])
{
  __m128d r;
  cell_wrap_3 c[6];
  cell_wrap_3 uA_data[6];
  cell_wrap_3 uB_data[6];
  double dIdx_data[6];
  double icA[6];
  double locbOfUA_data[6];
  double y_data[6];
  double b_y_data[5];
  double b_dv[2];
  double b_d;
  int idx[6];
  int tmp_data[6];
  int b_obj_Value_size[2];
  int obj_Value_size[2];
  int b_i;
  int b_idx;
  int b_loop_ub;
  int i;
  int j;
  int loop_ub;
  int nz;
  int tmp_size;
  int vectorUB;
  unsigned char b2a[7];
  signed char ii_data[6];
  boolean_T d[6];
  boolean_T uAInUB_data[6];
  boolean_T exitg1;
  boolean_T y;
  introsort(bnames, idx);
  c[0].f1.size[0] = 1;
  loop_ub = bnames[idx[0] - 1].f1.size[1];
  c[0].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[0].f1.data[i] = bnames[idx[0] - 1].f1.data[i];
  }
  c[1].f1.size[0] = 1;
  loop_ub = bnames[idx[1] - 1].f1.size[1];
  c[1].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[1].f1.data[i] = bnames[idx[1] - 1].f1.data[i];
  }
  c[2].f1.size[0] = 1;
  loop_ub = bnames[idx[2] - 1].f1.size[1];
  c[2].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[2].f1.data[i] = bnames[idx[2] - 1].f1.data[i];
  }
  c[3].f1.size[0] = 1;
  loop_ub = bnames[idx[3] - 1].f1.size[1];
  c[3].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[3].f1.data[i] = bnames[idx[3] - 1].f1.data[i];
  }
  c[4].f1.size[0] = 1;
  loop_ub = bnames[idx[4] - 1].f1.size[1];
  c[4].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[4].f1.data[i] = bnames[idx[4] - 1].f1.data[i];
  }
  c[5].f1.size[0] = 1;
  loop_ub = bnames[idx[5] - 1].f1.size[1];
  c[5].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[5].f1.data[i] = bnames[idx[5] - 1].f1.data[i];
  }
  d[0] = true;
  nz = 1;
  for (i = 0; i < 5; i++) {
    y = !b_strcmp(c[i].f1.data, c[i].f1.size, c[i + 1].f1.data,
                  c[i + 1].f1.size);
    d[i + 1] = y;
    nz += y;
  }
  tmp_size = eml_find(d, tmp_data);
  for (i = 0; i < tmp_size; i++) {
    dIdx_data[i] = tmp_data[i];
  }
  for (b_i = 0; b_i < tmp_size; b_i++) {
    if (b_i + 1 != tmp_size) {
      b_d = dIdx_data[b_i + 1] - 1.0;
      loop_ub = (int)dIdx_data[b_i];
      if (b_d < loop_ub) {
        b_idx = 0;
      } else {
        b_loop_ub = (int)b_d - loop_ub;
        b_idx = b_loop_ub + 1;
        j = ((b_loop_ub + 1) / 2) << 1;
        vectorUB = j - 2;
        for (i = 0; i <= vectorUB; i += 2) {
          b_dv[0] = i;
          b_dv[1] = i + 1;
          r = _mm_loadu_pd(&b_dv[0]);
          _mm_storeu_pd(&b_y_data[i], _mm_add_pd(_mm_set1_pd(loop_ub), r));
        }
        for (i = j; i <= b_loop_ub; i++) {
          b_y_data[i] = (double)loop_ub + (double)i;
        }
      }
      loop_ub = idx[(int)b_y_data[0] - 1];
      for (i = 2; i <= b_idx; i++) {
        b_loop_ub = idx[(int)b_y_data[i - 1] - 1];
        if (loop_ub > b_loop_ub) {
          loop_ub = b_loop_ub;
        }
      }
      locbOfUA_data[b_i] = loop_ub;
    } else {
      loop_ub = (int)dIdx_data[b_i];
      b_idx = 7 - loop_ub;
      b_loop_ub = 6 - loop_ub;
      j = ((7 - loop_ub) / 2) << 1;
      vectorUB = j - 2;
      for (i = 0; i <= vectorUB; i += 2) {
        b_dv[0] = i;
        b_dv[1] = i + 1;
        r = _mm_loadu_pd(&b_dv[0]);
        _mm_storeu_pd(&y_data[i], _mm_add_pd(_mm_set1_pd(loop_ub), r));
      }
      for (i = j; i <= b_loop_ub; i++) {
        y_data[i] = (double)loop_ub + (double)i;
      }
      loop_ub = idx[(int)y_data[0] - 1];
      for (i = 2; i <= b_idx; i++) {
        b_loop_ub = idx[(int)y_data[i - 1] - 1];
        if (loop_ub > b_loop_ub) {
          loop_ub = b_loop_ub;
        }
      }
      locbOfUA_data[b_i] = loop_ub;
    }
  }
  for (i = 0; i < 6; i++) {
    y_data[i] = d[i];
  }
  for (i = 0; i < 5; i++) {
    y_data[i + 1] += y_data[i];
  }
  for (i = 0; i < 6; i++) {
    icA[i] = y_data[i];
  }
  for (i = 0; i < 6; i++) {
    y_data[idx[i] - 1] = icA[i];
  }
  for (i = 0; i < tmp_size; i++) {
    uA_data[i].f1.size[0] = 1;
    loop_ub = (int)locbOfUA_data[i];
    b_loop_ub = bnames[loop_ub - 1].f1.size[1];
    uA_data[i].f1.size[1] = b_loop_ub;
    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      uA_data[i].f1.data[b_i] = bnames[loop_ub - 1].f1.data[b_i];
    }
  }
  b_loop_ub = cellstr_unique(anames, uB_data, dIdx_data, &loop_ub);
  memset(&uAInUB_data[0], 0, (unsigned int)nz * sizeof(boolean_T));
  memset(&locbOfUA_data[0], 0, (unsigned int)nz * sizeof(double));
  if (b_loop_ub > 0) {
    j = 0;
    vectorUB = 0;
    exitg1 = false;
    while ((!exitg1) && (vectorUB <= nz - 1)) {
      char b_obj_Value_data[18];
      char obj_Value_data[18];
      obj_Value_size[0] = 1;
      b_idx = uA_data[vectorUB].f1.size[1];
      obj_Value_size[1] = b_idx;
      for (i = 0; i < b_idx; i++) {
        obj_Value_data[i] = uA_data[vectorUB].f1.data[i];
      }
      b_obj_Value_size[0] = 1;
      loop_ub = uB_data[j].f1.size[1];
      b_obj_Value_size[1] = loop_ub;
      for (i = 0; i < loop_ub; i++) {
        b_obj_Value_data[i] = uB_data[j].f1.data[i];
      }
      if (!string_lt(obj_Value_data, obj_Value_size, b_obj_Value_data,
                     b_obj_Value_size)) {
        if (b_strcmp(uA_data[vectorUB].f1.data, uA_data[vectorUB].f1.size,
                     uB_data[j].f1.data, uB_data[j].f1.size)) {
          uAInUB_data[vectorUB] = true;
          locbOfUA_data[vectorUB] = dIdx_data[j];
          j++;
        } else {
          boolean_T exitg2;
          exitg2 = false;
          while ((!exitg2) && (j + 1 <= b_loop_ub)) {
            obj_Value_size[0] = 1;
            for (i = 0; i < b_idx; i++) {
              obj_Value_data[i] = uA_data[vectorUB].f1.data[i];
            }
            b_obj_Value_size[0] = 1;
            loop_ub = uB_data[j].f1.size[1];
            b_obj_Value_size[1] = loop_ub;
            for (i = 0; i < loop_ub; i++) {
              b_obj_Value_data[i] = uB_data[j].f1.data[i];
            }
            if (string_gt(obj_Value_data, obj_Value_size, b_obj_Value_data,
                          b_obj_Value_size)) {
              j++;
            } else {
              exitg2 = true;
            }
          }
          if ((j + 1 <= b_loop_ub) &&
              b_strcmp(uA_data[vectorUB].f1.data, uA_data[vectorUB].f1.size,
                       uB_data[j].f1.data, uB_data[j].f1.size)) {
            uAInUB_data[vectorUB] = true;
            locbOfUA_data[vectorUB] = dIdx_data[j];
            j++;
          }
        }
      }
      if (j + 1 > b_loop_ub) {
        exitg1 = true;
      } else {
        vectorUB++;
      }
    }
  }
  for (i = 0; i < 7; i++) {
    b2a[i] = 0U;
  }
  for (i = 0; i < 6; i++) {
    unsigned char u;
    loop_ub = (int)y_data[i] - 1;
    b_d = rt_roundd_snf(locbOfUA_data[loop_ub]);
    if (b_d < 256.0) {
      if (b_d >= 0.0) {
        u = (unsigned char)b_d;
      } else {
        u = 0U;
      }
    } else if (b_d >= 256.0) {
      u = MAX_uint8_T;
    } else {
      u = 0U;
    }
    b2a[i + 1] = u;
    d[i] = !uAInUB_data[loop_ub];
  }
  b_idx = 0;
  loop_ub = 0;
  exitg1 = false;
  while ((!exitg1) && (loop_ub < 6)) {
    if (d[loop_ub]) {
      b_idx++;
      ii_data[b_idx - 1] = (signed char)(loop_ub + 1);
      if (b_idx >= 6) {
        exitg1 = true;
      } else {
        loop_ub++;
      }
    } else {
      loop_ub++;
    }
  }
  if (b_idx < 1) {
    b_idx = 0;
  }
  y = true;
  loop_ub = 0;
  exitg1 = false;
  while ((!exitg1) && (loop_ub < 6)) {
    if (!uAInUB_data[(int)y_data[loop_ub] - 1]) {
      y = false;
      exitg1 = true;
    } else {
      loop_ub++;
    }
  }
  if (!y) {
    if (b_idx < 1) {
      b_idx = 0;
    } else {
      loop_ub = b_idx - 1;
      b_loop_ub = (b_idx / 2) << 1;
      j = b_loop_ub - 2;
      for (i = 0; i <= j; i += 2) {
        b_dv[0] = i;
        b_dv[1] = (double)i + 1.0;
        r = _mm_loadu_pd(&b_dv[0]);
        _mm_storeu_pd(&y_data[i], _mm_add_pd(_mm_set1_pd(1.0), r));
      }
      for (i = b_loop_ub; i <= loop_ub; i++) {
        y_data[i] = (double)i + 1.0;
      }
    }
    for (i = 0; i < b_idx; i++) {
      b2a[ii_data[i]] = (unsigned char)((unsigned char)y_data[i] + 6);
    }
  }
  for (i = 0; i < 6; i++) {
    bcodesout[i] = b2a[i + 1];
  }
}

/* End of code generation (convertCodes.c) */
