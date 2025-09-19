/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * introsort.c
 *
 * Code generation for function 'introsort'
 *
 */

/* Include files */
#include "introsort.h"
#include "insertionsort.h"
#include "predict_exercise_index_internal_types.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>

/* Function Definitions */
void introsort(const cell_wrap_3 cmp_workspace_c[6], int x[6])
{
  static const int offsets[4] = {0, 1, 2, 3};
  _mm_storeu_si128(
      (__m128i *)&x[0],
      _mm_add_epi32(
          _mm_set1_epi32(1),
          _mm_add_epi32(_mm_set1_epi32(0),
                        _mm_loadu_si128((const __m128i *)&offsets[0]))));
  x[4] = 5;
  x[5] = 6;
  insertionsort(x, cmp_workspace_c);
}

/* End of code generation (introsort.c) */
