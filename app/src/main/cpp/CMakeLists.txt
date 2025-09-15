# CMake의 최소 요구 버전을 설정합니다.
cmake_minimum_required(VERSION 3.22.1)

project("geunhwang")

# --- ▼▼▼ [추가] OpenMP 라이브러리를 찾고, 컴파일러에 활성화 옵션을 추가합니다. ▼▼▼ ---
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    message(STATUS "OpenMP found, enabling parallel support.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_library(
        geunhwang
        SHARED
        native-lib.cpp

        # (이전에 추가했던 모든 .cpp 파일 목록)
        abs.cpp
        blockedSummation.cpp
        categorical.cpp
        ClassificationKNN.cpp
        corrcoef.cpp
        eml_setop.cpp
        feature_extractor_codegen.cpp
        FFTImplementationCallback.cpp
        find.cpp
        findpeaks.cpp
        geomean.cpp
        gradient.cpp
        heapsort.cpp
        insertionsort.cpp
        introsort.cpp
        mean.cpp
        minOrMax.cpp
        norm.cpp
        pca.cpp
        pdist2.cpp
        periodogram.cpp
        predictExercise.cpp
        predictExercise_data.cpp
        predictExercise_initialize.cpp
        predictExercise_terminate.cpp
        rms.cpp
        rtGetInf.cpp
        rtGetNaN.cpp
        rt_nonfinite.cpp
        scanfornan.cpp
        sort.cpp
        sortIdx.cpp
        std.cpp
        strcmp.cpp
        strtrim.cpp
        sum.cpp
        topkrows.cpp
        trapz.cpp
        var.cpp
        xaxpy.cpp
        xdotc.cpp
        xnrm2.cpp
        xrot.cpp
        xrotg.cpp
        xswap.cpp
        xzlangeM.cpp
        xzlascl.cpp
        xzsvdc.cpp
)

find_library(
        log-lib
        log)

# --- ▼▼▼ [수정] 최종 라이브러리에 OpenMP와 log 라이브러리를 모두 연결합니다. ▼▼▼ ---
target_link_libraries(
        geunhwang
        OpenMP::OpenMP_CXX
        ${log-lib}
)