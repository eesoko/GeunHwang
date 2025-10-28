/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * updateCache.c
 *
 * Code generation for function 'updateCache'
 *
 */

/* Include files */
#include "updateCache.h"
#include "predict_exercise_index_data.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>

/* Function Definitions */
void updateCache(const double learnerscore[6], double cachedScore[6],
                 double *cachedWeights, boolean_T *cached,
                 const char combinerName[15], boolean_T obsIndices,
                 double score[6])
{
  int i;
  int i1;
  for (i = 0; i < 6; i++) {
    score[i] = cachedScore[i];
  }
  if (!*cached) {
    double cachedScore_data[6];
    double d;
    int trueCount;
    boolean_T b_bool;
    *cached = true;
    trueCount = 0;
    if (obsIndices) {
      trueCount = 1;
    }
    for (i = 0; i < 6; i++) {
      if (trueCount - 1 >= 0) {
        cachedScore_data[trueCount * i] = cachedScore[i] + learnerscore[i];
      }
    }
    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < trueCount; i1++) {
        cachedScore[i] = cachedScore_data[trueCount * i];
      }
    }
    d = *cachedWeights;
    if (trueCount - 1 >= 0) {
      d = *cachedWeights + 1.0;
    }
    *cachedWeights = d;
    b_bool = false;
    trueCount = 0;
    int exitg1;
    do {
      exitg1 = 0;
      if (trueCount < 15) {
        if (cv3[(int)combinerName[trueCount]] != cv3[(int)cv4[trueCount]]) {
          exitg1 = 1;
        } else {
          trueCount++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
    if (b_bool) {
      if (d == 0.0) {
        boolean_T exitg2;
        b_bool = true;
        trueCount = 0;
        exitg2 = false;
        while ((!exitg2) && (trueCount < 6)) {
          if (!(cachedScore[trueCount] == 0.0)) {
            b_bool = false;
            exitg2 = true;
          } else {
            trueCount++;
          }
        }
        if (b_bool) {
          for (i = 0; i < 6; i++) {
            score[i] = rtNaN;
          }
        } else {
          b_bool = true;
          trueCount = 0;
          exitg2 = false;
          while ((!exitg2) && (trueCount < 6)) {
            if (!(cachedScore[trueCount] < 0.0)) {
              b_bool = false;
              exitg2 = true;
            } else {
              trueCount++;
            }
          }
          if (b_bool) {
            for (i = 0; i < 6; i++) {
              score[i] = rtMinusInf;
            }
          } else {
            for (i = 0; i < 6; i++) {
              score[i] = rtInf;
            }
          }
        }
      } else {
        __m128d r;
        r = _mm_loadu_pd(&cachedScore[0]);
        _mm_storeu_pd(&score[0], _mm_div_pd(r, _mm_set1_pd(*cachedWeights)));
        r = _mm_loadu_pd(&cachedScore[2]);
        _mm_storeu_pd(&score[2], _mm_div_pd(r, _mm_set1_pd(*cachedWeights)));
        r = _mm_loadu_pd(&cachedScore[4]);
        _mm_storeu_pd(&score[4], _mm_div_pd(r, _mm_set1_pd(*cachedWeights)));
      }
    } else {
      for (i = 0; i < 6; i++) {
        score[i] = cachedScore[i];
      }
    }
  }
}

/* End of code generation (updateCache.c) */
