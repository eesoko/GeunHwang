// =========================================================================
// native-lib.cpp (최종 수정 완료)
// =========================================================================
#include <jni.h>
#include <string>
#include <vector>

// MATLAB에서 생성된 C 코드의 헤더 파일들을 포함합니다.
// 필요한 모든 타입 정의가 포함되도록 헤더를 추가했습니다.
extern "C" {
#include "feature_extractor_codegen.h"
#include "feature_extractor_codegen_types.h" // emx... 타입 정의 포함
#include "predict_exercise.h"
#include "predict_exercise_initialize.h"
#include "predict_exercise_terminate.h"
#include "predict_exercise_types.h"
}


extern "C" JNIEXPORT jstring JNICALL
// Kotlin의 private external 함수 이름에 맞춰 Native 접미사 추가
Java_com_example_geunhwang_presentation_ui_MainActivity_predictMotionNative(
        JNIEnv* env,
        jobject /* this */,
        jobjectArray sensorData, // 입력 1: 센서 데이터 (Nx6 크기의 2차원 배열)
        jdouble fs) {            // 입력 2: 샘플링 주파수 (double)

    // --- 1. Kotlin의 2차원 배열 -> C++ 벡터로 변환 ---
    int rows = env->GetArrayLength(sensorData);
    if (rows == 0) {
        return env->NewStringUTF("Input data is empty.");
    }
    jdoubleArray firstRow = (jdoubleArray)env->GetObjectArrayElement(sensorData, 0);
    int cols = env->GetArrayLength(firstRow);

    std::vector<double> raw_data_flat;
    raw_data_flat.reserve(rows * cols); // 메모리 미리 할당

    for (int i = 0; i < rows; i++) {
        jdoubleArray row = (jdoubleArray)env->GetObjectArrayElement(sensorData, i);
        double* row_elements = env->GetDoubleArrayElements(row, nullptr);
        raw_data_flat.insert(raw_data_flat.end(), row_elements, row_elements + cols);
        env->ReleaseDoubleArrayElements(row, row_elements, JNI_ABORT);
        env->DeleteLocalRef(row);
    }

    // --- 2. C 함수 호출을 위한 데이터 준비 ---
    int raw_data_size[2] = {rows, cols};
    double features[32]; // 특징 벡터는 고정 크기 배열
    categorical prediction_result; // 예측 결과는 categorical 구조체

    // --- 3. 생성된 C 함수들을 순서대로 호출 ---
    predict_exercise_initialize();

    // 3-1. 특징 추출기 호출 (emxArray 대신 C 기본 배열 타입 사용)
    feature_extractor_codegen(raw_data_flat.data(), raw_data_size, fs, features);

    // 3-2. 운동 예측 함수 호출
    predict_exercise(features, &prediction_result);

    predict_exercise_terminate();

    // --- 4. C의 categorical 구조체 결과 -> Kotlin의 String으로 변환 ---
    // 'codes'는 1부터 시작하는 인덱스이므로 -1을 해줍니다.
    int codeIndex = prediction_result.codes - 1;
    std::string result_str;
    if (codeIndex >= 0 && codeIndex < 6) { // 배열 범위 확인
        // categoryNames 배열에서 해당 인덱스의 문자열 데이터를 가져옵니다.
        cell_wrap_0 result_cell = prediction_result.categoryNames[codeIndex];
        result_str.assign(result_cell.f1.data, result_cell.f1.size[1]);
    } else {
        result_str = "Unknown";
    }

    // --- 5. 최종 결과 문자열 반환 ---
    // 예측된 운동 이름과 함께 임시로 횟수 "10"을 붙여서 반환합니다.
    // 추후 C++ 코드에서 횟수 계산 로직이 추가되면 이 부분을 수정해야 합니다.
    std::string final_output = result_str + ",10";

    return env->NewStringUTF(final_output.c_str());
}