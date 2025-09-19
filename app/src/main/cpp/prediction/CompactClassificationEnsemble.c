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
#include "predict_exercise_data.h"
#include "predict_exercise_internal_types.h"
#include "rt_nonfinite.h"
#include "strtrim.h"

/* Function Definitions */
int c_CompactClassificationEnsemble(const double Xin[32])
{
    // --- ▼▼▼ 이 부분이 핵심 수정 내용입니다. ▼▼▼
    double scores[6]; // 'unusedExpr' 대신 'scores'라는 이름으로 변경
    int i;
    int max_index = 0;
    double max_score = -1.0;

    // C 함수를 호출하여 6개 클래스에 대한 예측 점수를 얻습니다.
    CompactEnsemble_ensemblePredict(Xin, scores);

    // 가장 높은 점수를 가진 클래스의 인덱스를 찾습니다.
    for (i = 0; i < 6; i++) {
        if (scores[i] > max_score) {
            max_score = scores[i];
            max_index = i;
        }
    }

    // 가장 높은 점수의 인덱스를 반환합니다. (0~5 사이의 값)
    return max_index;
    // --- ▲▲▲ 여기까지가 핵심 수정 내용입니다. ▲▲▲
}

/* End of code generation (CompactClassificationEnsemble.c) */
