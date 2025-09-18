#include <jni.h>
#include <string>
#include <vector>

// MATLAB에서 생성된 C 코드의 헤더 파일들을 포함합니다.
// --- ▼▼▼ 모든 #include 경로에 "features/" 또는 "prediction/"을 추가했습니다. ▼▼▼
extern "C" {
#include "features/feature_extractor_codegen.h"
#include "features/feature_extractor_codegen_types.h"
#include "prediction/predict_exercise.h"
#include "prediction/predict_exercise_initialize.h"
#include "prediction/predict_exercise_terminate.h"
#include "prediction/predict_exercise_types.h"
#include "prediction/predict_exercise_internal_types.h"
#include "prediction/CompactClassificationModel.h"
}
// --- ▲▲▲ ---

typedef struct {
    double codes;
    cell_wrap_3 categoryNames[6];
} categorical;

extern "C" JNIEXPORT jstring JNICALL
Java_com_example_geunhwang_presentation_ui_MainActivity_predictMotionNative(
        JNIEnv* env,
        jobject /* this */,
        jobjectArray sensorData,
        jdouble fs) {

    // --- 1. Kotlin의 2차원 배열 -> C++ 벡터로 변환 ---
    int rows = env->GetArrayLength(sensorData);
    if (rows == 0) {
        return env->NewStringUTF("Input data is empty.");
    }
    jdoubleArray firstRow = (jdoubleArray)env->GetObjectArrayElement(sensorData, 0);
    int cols = env->GetArrayLength(firstRow);

    std::vector<double> raw_data_flat;
    raw_data_flat.reserve(rows * cols);

    for (int i = 0; i < rows; i++) {
        jdoubleArray row = (jdoubleArray)env->GetObjectArrayElement(sensorData, i);
        double* row_elements = env->GetDoubleArrayElements(row, nullptr);
        raw_data_flat.insert(raw_data_flat.end(), row_elements, row_elements + cols);
        env->ReleaseDoubleArrayElements(row, row_elements, JNI_ABORT);
        env->DeleteLocalRef(row);
    }

    // --- 2. C 함수 호출을 위한 데이터 준비 ---
    int raw_data_size[2] = {rows, cols};
    double features[32];
    categorical prediction_result;

    // --- 3. 생성된 C 함수들을 순서대로 호출 ---
    predict_exercise_initialize();
    feature_extractor_codegen(raw_data_flat.data(), raw_data_size, fs, features);
    predict_exercise(features, &prediction_result);
    predict_exercise_terminate();

    // --- 4. C의 categorical 구조체 결과 -> Kotlin의 String으로 변환 ---
    int codeIndex = prediction_result.codes - 1; // 'codes'는 1부터 시작
    std::string result_str;
    if (codeIndex >= 0 && codeIndex < 6) { // 배열 범위 확인
        cell_wrap_3 result_cell = prediction_result.categoryNames[codeIndex];
        result_str.assign(result_cell.f1.data, result_cell.f1.size[1]);
    } else {
        result_str = "Unknown";
    }

    // --- 5. 최종 결과 문자열 반환 ---
    std::string final_output = result_str + ",10"; // 임시로 횟수 "10" 추가

    return env->NewStringUTF(final_output.c_str());
}