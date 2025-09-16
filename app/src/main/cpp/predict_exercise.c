/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * predict_exercise.c
 *
 * Code generation for function 'predict_exercise'
 *
 */

/* Include files */
#include "predict_exercise.h"
#include "ClassificationECOC.h"
#include "predict_exercise_data.h"
#include "predict_exercise_initialize.h"
#include "predict_exercise_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void predict_exercise(const double input_features[32], categorical *label)
{
  if (!isInitialized_predict_exercise) {
    predict_exercise_initialize();
  }
  /*  --- predict_exercise.m 파일의 최종 내용 --- */
  /*  '깨끗한' 모델 파일인 cleanModel.mat을 로드합니다. */
  /*  이제 s는 바로 'actual_model' 객체를 포함하고 있습니다. */
  label->codes =
      ClassificationECOC_predict(input_features, label->categoryNames);
}

/* End of code generation (predict_exercise.c) */
