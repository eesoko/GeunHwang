#include <jni.h>
#include <string>
#include <vector>
#include <android/log.h>

// MATLAB Coder가 생성한 헤더 파일들을 포함합니다.
#include "predictExercise.h"
#include "predictExercise_types.h"
#include "rtwtypes.h"
#include "categorical.h"

extern "C" {

JNIEXPORT jstring JNICALL
Java_com_example_geunhwang_presentation_ui_ExerciseClassifier_classify(
        JNIEnv* env,
        jobject /* this */,
        jobjectArray sensorData) {

    __android_log_print(ANDROID_LOG_DEBUG, "GeunhwangCPP", "Real classification started!");

    // 1. Kotlin의 2D Array를 C++의 2D Vector로 변환합니다.
    int rows = env->GetArrayLength(sensorData);
    if (rows == 0) {
        return env->NewStringUTF("Error:EmptyData");
    }
    jfloatArray firstRow = (jfloatArray)env->GetObjectArrayElement(sensorData, 0);
    int cols = env->GetArrayLength(firstRow);
    env->DeleteLocalRef(firstRow);

    std::vector<std::vector<double>> sensorDataVector(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        jfloatArray row = (jfloatArray)env->GetObjectArrayElement(sensorData, i);
        jfloat* rowElements = env->GetFloatArrayElements(row, nullptr);
        for (int j = 0; j < cols; ++j) {
            sensorDataVector[i][j] = static_cast<double>(rowElements[j]);
        }
        env->ReleaseFloatArrayElements(row, rowElements, 0);
        env->DeleteLocalRef(row);
    }

    // 2. C++ Vector를 MATLAB Coder가 이해하는 coder::array 형태로 변환합니다.
    coder::array<double, 2U> matlabSensorData;
    matlabSensorData.set_size(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matlabSensorData[i + rows * j] = sensorDataVector[i][j];
        }
    }

    // 3. ✨✨✨ 실제 MATLAB AI 모델 함수를 호출합니다! ✨✨✨
    coder::categorical predictedExercise_cat;
    predictExercise(matlabSensorData, &predictedExercise_cat);

    // 4. ▼▼▼ [최종 수정] categorical 객체의 정확한 내부 구조에 접근합니다. ▼▼▼
    std::string predictedExercise_str;

    unsigned char predicted_code = predictedExercise_cat.codes;

    if (predicted_code > 0 && predicted_code <= 6) {
        coder::bounded_array<char, 18U, 2U> exerciseNameChars;
        exerciseNameChars = predictedExercise_cat.categoryNames[predicted_code - 1].f1;

        // 4.1. [수정] .size는 함수가 아닌 배열 변수이므로, 두 번째 요소(문자열 길이)에 접근합니다.
        int string_length = exerciseNameChars.size[1];

        // 4.2. [수정] [] 대신 .data 내부 배열에 접근합니다.
        for (int i = 0; i < string_length; ++i) {
            predictedExercise_str += exerciseNameChars.data[i];
        }
    } else {
        predictedExercise_str = "Unknown";
    }

    // 5. 최종 결과를 "Exercise:결과,Reps:12" 형태로 만들어 Kotlin에 반환합니다.
    std::string result = "Exercise:" + predictedExercise_str + ",Reps:12";

    __android_log_print(ANDROID_LOG_DEBUG, "GeunhwangCPP", "Predicted: %s", predictedExercise_str.c_str());

    return env->NewStringUTF(result.c_str());
}

} // extern "C"