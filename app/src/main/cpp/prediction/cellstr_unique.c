/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * cellstr_unique.c
 *
 * Code generation for function 'cellstr_unique'
 *
 */

/* Include files */
#include "cellstr_unique.h"
#include "find.h"
#include "introsort.h"
#include "predict_exercise_index_internal_types.h"
#include "rt_nonfinite.h"
#include "strcmp.h"
#include <emmintrin.h>

/* Function Definitions */
int b_cellstr_unique(const cell_wrap_3 a[6], cell_wrap_3 u_data[])
{
  cell_wrap_3 c[6];
  double y_data[6];
  double b_y_data[5];
  double b_dv[2];
  int dIdx_data[6];
  int ia_data[6];
  int idx[6];
  int b_i;
  int dIdx_size;
  int i;
  int loop_ub;
  int u_size;
  boolean_T d[6];
  introsort(a, idx);
  c[0].f1.size[0] = 1;
  loop_ub = a[idx[0] - 1].f1.size[1];
  c[0].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[0].f1.data[i] = a[idx[0] - 1].f1.data[i];
  }
  c[1].f1.size[0] = 1;
  loop_ub = a[idx[1] - 1].f1.size[1];
  c[1].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[1].f1.data[i] = a[idx[1] - 1].f1.data[i];
  }
  c[2].f1.size[0] = 1;
  loop_ub = a[idx[2] - 1].f1.size[1];
  c[2].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[2].f1.data[i] = a[idx[2] - 1].f1.data[i];
  }
  c[3].f1.size[0] = 1;
  loop_ub = a[idx[3] - 1].f1.size[1];
  c[3].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[3].f1.data[i] = a[idx[3] - 1].f1.data[i];
  }
  c[4].f1.size[0] = 1;
  loop_ub = a[idx[4] - 1].f1.size[1];
  c[4].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[4].f1.data[i] = a[idx[4] - 1].f1.data[i];
  }
  c[5].f1.size[0] = 1;
  loop_ub = a[idx[5] - 1].f1.size[1];
  c[5].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[5].f1.data[i] = a[idx[5] - 1].f1.data[i];
  }
  d[0] = true;
  u_size = 1;
  for (i = 0; i < 5; i++) {
    boolean_T b;
    b = !b_strcmp(c[i].f1.data, c[i].f1.size, c[i + 1].f1.data,
                  c[i + 1].f1.size);
    d[i + 1] = b;
    u_size += b;
  }
  dIdx_size = eml_find(d, dIdx_data);
  for (b_i = 0; b_i < dIdx_size; b_i++) {
    int minval;
    if (b_i + 1 != dIdx_size) {
      double b_d;
      int y_size_idx_1;
      loop_ub = dIdx_data[b_i];
      b_d = (double)dIdx_data[b_i + 1] - 1.0;
      if (b_d < loop_ub) {
        y_size_idx_1 = 0;
      } else {
        int scalarLB;
        int vectorUB;
        minval = (int)b_d - loop_ub;
        y_size_idx_1 = minval + 1;
        scalarLB = ((minval + 1) / 2) << 1;
        vectorUB = scalarLB - 2;
        for (i = 0; i <= vectorUB; i += 2) {
          __m128d r;
          b_dv[0] = i;
          b_dv[1] = i + 1;
          r = _mm_loadu_pd(&b_dv[0]);
          _mm_storeu_pd(&b_y_data[i], _mm_add_pd(_mm_set1_pd(loop_ub), r));
        }
        for (i = scalarLB; i <= minval; i++) {
          b_y_data[i] = (double)loop_ub + (double)i;
        }
      }
      minval = idx[(int)b_y_data[0] - 1];
      for (i = 2; i <= y_size_idx_1; i++) {
        loop_ub = idx[(int)b_y_data[i - 1] - 1];
        if (minval > loop_ub) {
          minval = loop_ub;
        }
      }
      ia_data[b_i] = minval;
    } else {
      int scalarLB;
      int vectorUB;
      int y_size_idx_1;
      loop_ub = dIdx_data[b_i];
      y_size_idx_1 = 7 - loop_ub;
      minval = 6 - loop_ub;
      scalarLB = ((7 - loop_ub) / 2) << 1;
      vectorUB = scalarLB - 2;
      for (i = 0; i <= vectorUB; i += 2) {
        __m128d r;
        b_dv[0] = i;
        b_dv[1] = i + 1;
        r = _mm_loadu_pd(&b_dv[0]);
        _mm_storeu_pd(&y_data[i], _mm_add_pd(_mm_set1_pd(loop_ub), r));
      }
      for (i = scalarLB; i <= minval; i++) {
        y_data[i] = (double)loop_ub + (double)i;
      }
      minval = idx[(int)y_data[0] - 1];
      for (i = 2; i <= y_size_idx_1; i++) {
        loop_ub = idx[(int)y_data[i - 1] - 1];
        if (minval > loop_ub) {
          minval = loop_ub;
        }
      }
      ia_data[b_i] = minval;
    }
    u_data[b_i].f1.size[0] = 1;
    loop_ub = a[minval - 1].f1.size[1];
    u_data[b_i].f1.size[1] = loop_ub;
    for (i = 0; i < loop_ub; i++) {
      u_data[b_i].f1.data[i] = a[ia_data[b_i] - 1].f1.data[i];
    }
  }
  return u_size;
}

int cellstr_unique(const cell_wrap_3 a[6], cell_wrap_3 u_data[],
                   double ia_data[], int *ia_size)
{
  cell_wrap_3 c[6];
  double y_data[6];
  double b_y_data[5];
  double b_dv[2];
  int dIdx_data[6];
  int idx[6];
  int b_i;
  int dIdx_size;
  int i;
  int loop_ub;
  int u_size;
  boolean_T d[6];
  introsort(a, idx);
  c[0].f1.size[0] = 1;
  loop_ub = a[idx[0] - 1].f1.size[1];
  c[0].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[0].f1.data[i] = a[idx[0] - 1].f1.data[i];
  }
  c[1].f1.size[0] = 1;
  loop_ub = a[idx[1] - 1].f1.size[1];
  c[1].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[1].f1.data[i] = a[idx[1] - 1].f1.data[i];
  }
  c[2].f1.size[0] = 1;
  loop_ub = a[idx[2] - 1].f1.size[1];
  c[2].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[2].f1.data[i] = a[idx[2] - 1].f1.data[i];
  }
  c[3].f1.size[0] = 1;
  loop_ub = a[idx[3] - 1].f1.size[1];
  c[3].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[3].f1.data[i] = a[idx[3] - 1].f1.data[i];
  }
  c[4].f1.size[0] = 1;
  loop_ub = a[idx[4] - 1].f1.size[1];
  c[4].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[4].f1.data[i] = a[idx[4] - 1].f1.data[i];
  }
  c[5].f1.size[0] = 1;
  loop_ub = a[idx[5] - 1].f1.size[1];
  c[5].f1.size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    c[5].f1.data[i] = a[idx[5] - 1].f1.data[i];
  }
  d[0] = true;
  u_size = 1;
  for (i = 0; i < 5; i++) {
    boolean_T b;
    b = !b_strcmp(c[i].f1.data, c[i].f1.size, c[i + 1].f1.data,
                  c[i + 1].f1.size);
    d[i + 1] = b;
    u_size += b;
  }
  dIdx_size = eml_find(d, dIdx_data);
  *ia_size = dIdx_size;
  for (b_i = 0; b_i < dIdx_size; b_i++) {
    int minval;
    if (b_i + 1 != dIdx_size) {
      double b_d;
      int y_size_idx_1;
      loop_ub = dIdx_data[b_i];
      b_d = (double)dIdx_data[b_i + 1] - 1.0;
      if (b_d < loop_ub) {
        y_size_idx_1 = 0;
      } else {
        int scalarLB;
        int vectorUB;
        minval = (int)b_d - loop_ub;
        y_size_idx_1 = minval + 1;
        scalarLB = ((minval + 1) / 2) << 1;
        vectorUB = scalarLB - 2;
        for (i = 0; i <= vectorUB; i += 2) {
          __m128d r;
          b_dv[0] = i;
          b_dv[1] = i + 1;
          r = _mm_loadu_pd(&b_dv[0]);
          _mm_storeu_pd(&b_y_data[i], _mm_add_pd(_mm_set1_pd(loop_ub), r));
        }
        for (i = scalarLB; i <= minval; i++) {
          b_y_data[i] = (double)loop_ub + (double)i;
        }
      }
      minval = idx[(int)b_y_data[0] - 1];
      for (i = 2; i <= y_size_idx_1; i++) {
        loop_ub = idx[(int)b_y_data[i - 1] - 1];
        if (minval > loop_ub) {
          minval = loop_ub;
        }
      }
      ia_data[b_i] = minval;
    } else {
      int scalarLB;
      int vectorUB;
      int y_size_idx_1;
      loop_ub = dIdx_data[b_i];
      y_size_idx_1 = 7 - loop_ub;
      minval = 6 - loop_ub;
      scalarLB = ((7 - loop_ub) / 2) << 1;
      vectorUB = scalarLB - 2;
      for (i = 0; i <= vectorUB; i += 2) {
        __m128d r;
        b_dv[0] = i;
        b_dv[1] = i + 1;
        r = _mm_loadu_pd(&b_dv[0]);
        _mm_storeu_pd(&y_data[i], _mm_add_pd(_mm_set1_pd(loop_ub), r));
      }
      for (i = scalarLB; i <= minval; i++) {
        y_data[i] = (double)loop_ub + (double)i;
      }
      minval = idx[(int)y_data[0] - 1];
      for (i = 2; i <= y_size_idx_1; i++) {
        loop_ub = idx[(int)y_data[i - 1] - 1];
        if (minval > loop_ub) {
          minval = loop_ub;
        }
      }
      ia_data[b_i] = minval;
    }
    u_data[b_i].f1.size[0] = 1;
    loop_ub = a[minval - 1].f1.size[1];
    u_data[b_i].f1.size[1] = loop_ub;
    for (i = 0; i < loop_ub; i++) {
      u_data[b_i].f1.data[i] = a[(int)ia_data[b_i] - 1].f1.data[i];
    }
  }
  return u_size;
}

/* End of code generation (cellstr_unique.c) */
