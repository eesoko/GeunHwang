#include <jni.h>
#include <string>
#include <vector>

// MATLAB에서 생성된 C 코드의 헤더 파일들을 포함합니다.
extern "C" {
#include "features/feature_extractor_codegen.h"
#include "prediction/predict_exercise.h"
#include "prediction/predict_exercise_initialize.h"
#include "prediction/predict_exercise_terminate.h"
}


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
    double features[32];

    // --- 3. 생성된 C 함수들을 순서대로 호출 ---
    predict_exercise_initialize();
    feature_extractor_codegen(raw_data_flat.data(), new int[2]{rows, cols}, fs, features);

    // predict_exercise()를 호출하고 예측된 운동의 인덱스(0~5)를 받습니다.
    int predicted_index = predict_exercise(features);

    predict_exercise_terminate();

    // --- 4. 반환된 인덱스 -> 실제 운동 이름(String)으로 변환 ---
    std::string result_str;
    // MATLAB에서 정의한 운동 이름 순서대로 배열을 만듭니다.
    const char* exercise_names[6] = {
            "Dumbbell Curl",
            "Lunge",
            "Overhead Press",
            "Push Up",
            "Side Lateral Raise",
            "Squat"
    };

    if (predicted_index >= 0 && predicted_index < 6) {
        result_str = exercise_names[predicted_index];
    } else {
        result_str = "Unknown"; // 예외 처리
    }

    // --- 5. 최종 결과 문자열 반환 ---
    // 예시: "Squat,10"
    std::string final_output = result_str + ",10"; // 임시로 횟수 "10" 추가

    return env->NewStringUTF(final_output.c_str());
}