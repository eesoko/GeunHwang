cmake_minimum_required(VERSION 3.4.1)

# 1. 프로젝트에서 C++17 표준을 사용하도록 설정합니다.
set(CMAKE_CXX_STANDARD 17)

# 2. OpenMP 라이브러리를 찾고, 컴파일러에 OpenMP 지원을 활성화합니다.
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    message(STATUS "OpenMP found, enabling flags.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
else()
    message(WARNING "OpenMP not found.")
endif()

add_library(
    native-lib
    SHARED
    native-lib.cpp
    abs.c
    blockedSummation.c
    ClassificationECOC.c
    CompactClassificationModel.c
    corrcoef.c
    eml_setop.c
    feature_extractor_codegen.c
    feature_extractor_codegen_data.c
    feature_extractor_codegen_emxutil.c
    feature_extractor_codegen_initialize.c
    feature_extractor_codegen_terminate.c
    FFTImplementationCallback.c
    findpeaks.c
    geomean.c
    gradient.c
    heapsort.c
    insertionsort.c
    introsort.c
    loss.c
    mean.c
    minOrMax.c
    norm.c
    pca.c
    periodogram.c
    Poly.c
    predict_exercise.c
    predict_exercise_data.c
    predict_exercise_initialize.c
    predict_exercise_terminate.c
    rms.c
    rtGetInf.c
    rtGetNaN.c
    rt_nonfinite.c
    sort.c
    sortIdx.c
    sqrt.c
    std.c
    strtrim.c
    sum.c
    trapz.c
    tmwtypes.h
    var.c
    xaxpy.c
    xdotc.c
    xnrm2.c
    xrot.c
    xrotg.c
    xswap.c
    xzlangeM.c
    xzlascl.c
    xzsvdc.c
)

find_library(
    log-lib
    log
)

target_link_libraries(
    native-lib
    ${log-lib}

    OpenMP::OpenMP_C
    OpenMP::OpenMP_CXX
)
