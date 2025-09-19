cmake_minimum_required(VERSION 3.10)
project(GeunHwangNative LANGUAGES C CXX)

add_library(native-lib SHARED
    native-lib.cpp
    features/FFTImplementationCallback.c
    features/abs.c
    features/blockedSummation.c
    features/corrcoef.c
    features/eml_setop.c
    features/feature_extractor_codegen.c
    features/feature_extractor_codegen_data.c
    features/feature_extractor_codegen_emxAPI.c
    features/feature_extractor_codegen_emxutil.c
    features/feature_extractor_codegen_initialize.c
    features/feature_extractor_codegen_terminate.c
    features/findpeaks.c
    features/geomean.c
    features/gradient.c
    features/heapsort.c
    features/insertionsort.c
    features/introsort.c
    features/mean.c
    features/minOrMax.c
    features/norm.c
    features/pca.c
    features/periodogram.c
    features/rms.c
    features/rtGetInf.c
    features/rtGetNaN.c
    features/rt_nonfinite.c
    features/sort.c
    features/sortIdx.c
    features/sqrt.c
    features/std.c
    features/sum.c
    features/trapz.c
    features/var.c
    features/xaxpy.c
    features/xdotc.c
    features/xnrm2.c
    features/xrot.c
    features/xrotg.c
    features/xswap.c
    features/xzlangeM.c
    features/xzlascl.c
    features/xzsvdc.c
    prediction/CholeskyDiscriminant.c
    prediction/CompactClassificationDiscriminant.c
    prediction/CompactClassificationEnsemble.c
    prediction/CompactClassificationModel.c
    prediction/CompactEnsemble.c
    prediction/cellstr_unique.c
    prediction/convertCodes.c
    prediction/find.c
    prediction/ismember.c
    prediction/predictOneWithCache.c
    prediction/predict_exercise_index.c
    prediction/predict_exercise_index_data.c
    prediction/predict_exercise_index_initialize.c
    prediction/predict_exercise_index_terminate.c
    prediction/strcmp.c
    prediction/string1.c
    prediction/strtrim.c
    prediction/updateCache.c
)

find_library(log-lib log)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(native-lib ${log-lib} OpenMP::OpenMP_CXX)
else()
    message(FATAL_ERROR "OpenMP not found")
endif()

target_include_directories(native-lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
