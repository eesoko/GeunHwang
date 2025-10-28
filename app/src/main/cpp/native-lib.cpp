//
// Created by SeokHoLee on 2025-10-28.
//
// In: main/cpp/native-lib.cpp

#include <jni.h>
#include <string>
#include <vector>
#include <sstream>
#include <android/log.h>
#include <cmath>

// ... (include 문, getExerciseName 함수 등 상단 내용은 동일) ...
#include "features/rtwtypes.h"
#include "features/rt_nonfinite.h"
#include "features/feature_extractor_codegen_emxAPI.h"
#include "features/feature_extractor_codegen.h"
#include "features/feature_extractor_codegen_initialize.h"
#include "features/feature_extractor_codegen_terminate.h"
#include "prediction/predict_exercise_index.h"
#include "prediction/predict_exercise_index_initialize.h"
#include "prediction/predict_exercise_index_terminate.h"
#include "features/feature_extractor_codegen_types.h"

#define LOG_TAG "NativeLib"
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__)

// 운동 인덱스를 문자열로 변환하는 헬퍼 함수
std::string getExerciseName(int index) {
    switch (index) {
        case 1: return "Dumbbell Curl"; //
        case 2: return "Squat";         //
        case 3: return "Overhead Press"; //
        case 4: return "Push Up";        //
        case 5: return "Side Lateral Raise"; //
        case 6: return "Lunge";         //
        case 7: return "Dumbbell Row"; //
        default: return "Unknown";
    }
}

extern "C" JNIEXPORT jstring JNICALL
Java_com_example_geunhwang_presentation_ui_MainActivity_predictMotionNativeJNI(
        JNIEnv *env,
        jobject /* this */,
        jobjectArray sensorDataArray, // Kotlin의 List<FloatArray> -> Array<FloatArray>
        jdouble fs) {                 // Kotlin의 Double -> jdouble

    LOGD("Native function predictMotionNative called.");

    // 1. C++ 함수들 초기화
    feature_extractor_codegen_initialize();
    predict_exercise_index_initialize();

    // 2. JNI 데이터(jobjectArray)를 C++ emxArray로 변환
    int numRows = env->GetArrayLength(sensorDataArray);
    if (numRows == 0) {
        return env->NewStringUTF("Error: Empty data");
    }

    // emxArray 초기화 (N x 6)
    emxArray_real_T *sensorData = emxCreate_real_T(numRows, 6);

    // Kotlin의 List<FloatArray>를 C++의 emxArray로 복사
    for (int i = 0; i < numRows; ++i) {
        auto row = (jfloatArray) env->GetObjectArrayElement(sensorDataArray, i);
        jfloat *rowData = env->GetFloatArrayElements(row, nullptr);

        // 6축 데이터 (ax, ay, az, gx, gy, gz)
        sensorData->data[i] = rowData[0];           // ax
        sensorData->data[i + numRows] = rowData[1]; // ay
        sensorData->data[i + 2 * numRows] = rowData[2]; // az
        sensorData->data[i + 3 * numRows] = rowData[3]; // gx
        sensorData->data[i + 4 * numRows] = rowData[4]; // gy
        sensorData->data[i + 5 * numRows] = rowData[5]; // gz

        env->ReleaseFloatArrayElements(row, rowData, 0);
        env->DeleteLocalRef(row);
    }
    LOGD("Data converted. Rows: %d", numRows);

    // 3. 특징 추출기 실행
    // 32개 특징을 저장할 출력 emxArray 생성 (1 x 32)
    double features_array[32];

    // C++ 특징 추출기 함수 호출
    // void feature_extractor_codegen(const emxArray_real_T *raw_data, double Fs_actual, double features[32]);
    feature_extractor_codegen(sensorData, (double)fs, features_array); // emxArray 대신 배열 전달, real_T 대신 double 전달
    LOGD("Feature extraction complete. Feature[0]: %f", features_array[0]);

    // 4. 운동 분류기 실행
    int predicted_index_int; // 결과를 저장할 int 변수 (Header 파일에 맞춰 수정)

    // C++ 분류기 함수 호출
    // int predict_exercise_index(const double features[32]);
    predicted_index_int = predict_exercise_index(features_array); // 포인터 대신 반환값 받기, emxArray 대신 배열 전달
    LOGD("Prediction complete. Index: %d", predicted_index_int);

    // 5. 결과 조합 (운동 이름, 횟수)
    std::string exerciseName = getExerciseName(predicted_index_int);

    // 횟수(RPM)는 32개 특징 중 28번째 값(인덱스 27)입니다.
    double rpm = features_array[27];
    int reps = static_cast<int>(std::round(rpm));
    reps = std::max(1, reps); // 최소 1회 보장
    LOGD("Reps calculated: %d", reps);

    // 6. C++ 메모리 해제
    emxDestroyArray_real_T(sensorData);
    feature_extractor_codegen_terminate();
    predict_exercise_index_terminate();

    // 7. Kotlin으로 "운동이름,횟수" 포맷의 문자열 반환
    std::stringstream ss;
    ss << exerciseName << "," << reps;
    std::string resultString = ss.str();

    return env->NewStringUTF(resultString.c_str());
}