/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * predict_exercise_index.c
 *
 * Code generation for function 'predict_exercise_index'
 *
 */

/* Include files */
#include "predict_exercise_index.h"
#include "CompactClassificationEnsemble.h"
#include "minOrMax.h"
#include "predict_exercise_index_data.h"
#include "predict_exercise_index_initialize.h"
#include "rt_nonfinite.h"

/* Function Definitions */
int predict_exercise_index(const double features[32])
{
  double score[6];
  int idx;
  if (!isInitialized_predict_exercise_index) {
    predict_exercise_index_initialize();
  }
  /*  =========================================================================
   */
  /*  [ SCRIPT ]    : predict_exercise_index.m */
  /*  [ VERSION ]   : 1.0 */
  /*  [ PURPOSE ]   : Feature vector(1x32)를 입력받아 운동 라벨을 문자열로 반환
   */
  /*  =========================================================================
   */
  c_CompactClassificationEnsemble(features, score);
  /*  점수 행렬 */
  maximum(score, &idx);
  /*  가장 큰 확률 클래스 */
  /*  클래스 번호 (1, 2, 3, …) */
  return idx;
}

/* End of code generation (predict_exercise_index.c) */
