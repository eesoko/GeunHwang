# CMake 최소 요구 버전 설정
cmake_minimum_required(VERSION 3.22.1)

# 프로젝트 이름 설정 (경고 메시지 제거)
project("geunhwang_native")

# OpenMP 라이브러리 찾기 및 컴파일러 옵션 설정
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

# features 폴더와 prediction 폴더에 있는 모든 .c와 .cpp 파일을 자동으로 찾습니다.
file(GLOB_RECURSE SOURCES
        "features/*.c"
        "features/*.cpp"
        "prediction/*.c"
        "prediction/*.cpp"
)

# 네이티브 라이브러리를 빌드합니다.
add_library(
        native-lib
        SHARED

        # JNI 브릿지 파일
        native-lib.cpp

        # 위에서 자동으로 찾은 모든 소스 파일 목록
        ${SOURCES}
)

# --- ▼▼▼ 이 줄을 추가했습니다 ▼▼▼ ---
# 'native-lib'를 빌드할 때, 현재 CMakeLists.txt가 있는 폴더(즉, 'cpp' 폴더)를
# 헤더 검색 경로에 추가합니다.
target_include_directories(native-lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(native-lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/features)
target_include_directories(native-lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/prediction)
# --- ▲▲▲ ---

# 로그 라이브러리를 찾습니다.
find_library(
        log-lib
        log
)

# 라이브러리들을 최종 결과물에 연결(링크)합니다.
target_link_libraries(
        native-lib
        ${log-lib}

        # OpenMP 라이브러리 연결
        OpenMP::OpenMP_C
        OpenMP::OpenMP_CXX
)