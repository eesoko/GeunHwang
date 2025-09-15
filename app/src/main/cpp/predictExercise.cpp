//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// predictExercise.cpp
//
// Code generation for function 'predictExercise'
//

// Include files
#include "predictExercise.h"
#include "ClassificationKNN.h"
#include "abs.h"
#include "categorical.h"
#include "corrcoef.h"
#include "feature_extractor_codegen.h"
#include "findpeaks.h"
#include "geomean.h"
#include "gradient.h"
#include "mean.h"
#include "minOrMax.h"
#include "norm.h"
#include "pca.h"
#include "periodogram.h"
#include "predictExercise_data.h"
#include "predictExercise_initialize.h"
#include "rms.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "std.h"
#include "sum.h"
#include "trapz.h"
#include "var.h"
#include "coder_array.h"
#include <cmath>
#include <cstring>

// Variable Definitions
static bool modelObj_not_empty;

// Function Declarations
static void binary_expand_op_4(coder::array<double, 1U> &in1,
                               const coder::array<double, 1U> &in3,
                               const coder::array<double, 1U> &in4);

// Function Definitions
static void binary_expand_op_4(coder::array<double, 1U> &in1,
                               const coder::array<double, 1U> &in3,
                               const coder::array<double, 1U> &in4)
{
  coder::array<double, 1U> r;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  if (in4.size(0) == 1) {
    if (in1.size(0) == 1) {
      loop_ub = in3.size(0);
    } else {
      loop_ub = in1.size(0);
    }
  } else {
    loop_ub = in4.size(0);
  }
  r.set_size(loop_ub);
  stride_0_0 = (in3.size(0) != 1);
  stride_1_0 = (in1.size(0) != 1);
  stride_2_0 = (in4.size(0) != 1);
  for (int i{0}; i < loop_ub; i++) {
    double b_varargin_1;
    double c_varargin_1;
    double varargin_1;
    varargin_1 = in3[i * stride_0_0] * 50.0;
    b_varargin_1 = in1[i * stride_1_0] * 50.0;
    c_varargin_1 = in4[i * stride_2_0] * 50.0;
    r[i] = (varargin_1 * varargin_1 + b_varargin_1 * b_varargin_1) +
           c_varargin_1 * c_varargin_1;
  }
  in1.set_size(loop_ub);
  for (int i{0}; i < loop_ub; i++) {
    in1[i] = r[i];
  }
}

void predictExercise(const coder::array<double, 2U> &sensorData,
                     coder::categorical *predictedExercise)
{
  static coder::ClassificationKNN modelObj;
  coder::array<double, 2U> a;
  coder::array<double, 2U> a__5;
  coder::array<double, 2U> b_a;
  coder::array<double, 2U> r1;
  coder::array<double, 1U> a_perp;
  coder::array<double, 1U> accel_mag;
  coder::array<double, 1U> ax;
  coder::array<double, 1U> ay;
  coder::array<double, 1U> az;
  coder::array<double, 1U> gx;
  coder::array<double, 1U> gy;
  coder::array<double, 1U> gyro_mag;
  coder::array<double, 1U> gz;
  coder::array<double, 1U> jerk_mag;
  coder::array<double, 1U> r;
  coder::array<double, 1U> t;
  coder::array<int, 2U> y;
  coder::array<bool, 1U> b_gyro_mag;
  double features[32];
  double a__4_data[9];
  double dv[4];
  double g_vec[3];
  double latent_data[3];
  double b_varargin_1;
  double d;
  double d1;
  double d2;
  double g_norm;
  double varargin_1;
  int a__4_size[2];
  int loop_ub;
  int nx;
  if (!isInitialized_predictExercise) {
    predictExercise_initialize();
  }
  if (!modelObj_not_empty) {
    //  ▼▼▼ [수정] 새로운 Coder 전용 모델 로더와 파일 이름을 사용합니다. ▼▼▼
    modelObj.init();
    modelObj_not_empty = true;
  }
  //  =========================================================================
  //  [ SCRIPT ]    : feature_extractor_codegen.m
  //  [ VERSION ]   : 2.1 (Coder-Friendly Version)
  //  [ AUTHOR ]    : GeunHwang Project (Adapted by Gemini)
  //  [ DATE ]      : 2025-09-13
  //
  //  [ OVERVIEW ]
  //    'feature_extractor_elite.m'에서 C++ 코드 생성이 가능한 순수 계산
  //    로직만을 추출한 버전. 파일 입출력, UI, 지원되지 않는 함수 제거.
  //  =========================================================================
  //  --- 0. 입력 및 상수 정의 ---
  loop_ub = sensorData.size(0);
  ax.set_size(sensorData.size(0));
  ay.set_size(sensorData.size(0));
  az.set_size(sensorData.size(0));
  gx.set_size(sensorData.size(0));
  gy.set_size(sensorData.size(0));
  gz.set_size(sensorData.size(0));
  for (int k{0}; k < loop_ub; k++) {
    ax[k] = sensorData[k];
    ay[k] = sensorData[k + sensorData.size(0)];
    az[k] = sensorData[k + sensorData.size(0) * 2];
    gx[k] = sensorData[k + sensorData.size(0) * 3];
    gy[k] = sensorData[k + sensorData.size(0) * 4];
    gz[k] = sensorData[k + sensorData.size(0) * 5];
  }
  //  샘플링 주파수(가정)
  if (sensorData.size(0) - 1 < 0) {
    y.set_size(1, 0);
  } else {
    y.set_size(1, sensorData.size(0));
    for (int k{0}; k < loop_ub; k++) {
      y[k] = k;
    }
  }
  nx = y.size(1);
  t.set_size(y.size(1));
  for (int k{0}; k < nx; k++) {
    t[k] = static_cast<double>(y[k]) / 50.0;
  }
  //  출력 변수의 크기를 1x32로 미리 고정합니다.
  //  --- 기본 물리량 ---
  accel_mag.set_size(sensorData.size(0));
  gyro_mag.set_size(sensorData.size(0));
  for (int k{0}; k < loop_ub; k++) {
    g_norm = sensorData[k];
    varargin_1 = sensorData[k + sensorData.size(0)];
    b_varargin_1 = sensorData[k + sensorData.size(0) * 2];
    accel_mag[k] = (g_norm * g_norm + varargin_1 * varargin_1) +
                   b_varargin_1 * b_varargin_1;
    g_norm = sensorData[k + sensorData.size(0) * 3];
    varargin_1 = sensorData[k + sensorData.size(0) * 4];
    b_varargin_1 = sensorData[k + sensorData.size(0) * 5];
    gyro_mag[k] = (g_norm * g_norm + varargin_1 * varargin_1) +
                  b_varargin_1 * b_varargin_1;
  }
  nx = accel_mag.size(0);
  for (int k{0}; k < nx; k++) {
    accel_mag[k] = std::sqrt(accel_mag[k]);
  }
  nx = gyro_mag.size(0);
  for (int k{0}; k < nx; k++) {
    gyro_mag[k] = std::sqrt(gyro_mag[k]);
  }
  //  ==================== 기존 정예 특징 ====================
  //  기둥 1: 자세(Orientation)
  d = coder::mean(ax);
  features[0] = d;
  d1 = coder::mean(ay);
  features[1] = d1;
  d2 = coder::mean(az);
  features[2] = d2;
  varargin_1 = coder::b_std(ax);
  features[3] = varargin_1;
  b_varargin_1 = coder::b_std(ay);
  features[4] = b_varargin_1;
  features[5] = coder::b_std(az);
  //  기둥 2: 궤적(Trajectory)
  features[6] = coder::rms(gx);
  features[7] = coder::rms(gy);
  features[8] = coder::rms(gz);
  //  --- 유틸 함수 (코드 생성용) ---
  //  Coder가 지원하는 corrcoef를 사용하고, std가 0인 경우를 방어합니다.
  if ((varargin_1 == 0.0) || (b_varargin_1 == 0.0)) {
    features[9] = 0.0;
  } else {
    coder::corrcoef(ax, ay, dv);
    features[9] = dv[2];
  }
  //  --- 유틸 함수 (코드 생성용) ---
  //  Coder가 지원하는 corrcoef를 사용하고, std가 0인 경우를 방어합니다.
  if ((b_varargin_1 == 0.0) || (coder::b_std(gz) == 0.0)) {
    features[10] = 0.0;
  } else {
    coder::corrcoef(ay, gz, dv);
    features[10] = dv[2];
  }
  //  기둥 3: 주기/강도(Periodicity & Intensity)
  features[11] = coder::rms(accel_mag);
  if (accel_mag.size(0) == 0) {
    features[12] = 0.0;
  } else {
    features[12] = coder::internal::maximum(accel_mag) -
                   coder::internal::minimum(accel_mag);
  }
  features[13] = coder::rms(gyro_mag);
  coder::periodogram(accel_mag, jerk_mag, a_perp);
  coder::internal::maximum(jerk_mag, nx);
  features[14] = a_perp[nx - 1];
  coder::periodogram(gyro_mag, accel_mag, jerk_mag);
  coder::internal::maximum(accel_mag, nx);
  features[15] = jerk_mag[nx - 1];
  //  숨겨진 보석: 운동의 질(Quality)
  coder::gradient(ax, r);
  coder::gradient(ay, jerk_mag);
  coder::gradient(az, accel_mag);
  if (r.size(0) == 1) {
    nx = jerk_mag.size(0);
  } else {
    nx = r.size(0);
  }
  if ((r.size(0) == jerk_mag.size(0)) && (nx == accel_mag.size(0))) {
    nx = r.size(0);
    jerk_mag.set_size(r.size(0));
    for (int k{0}; k < nx; k++) {
      varargin_1 = r[k] * 50.0;
      b_varargin_1 = jerk_mag[k] * 50.0;
      g_norm = accel_mag[k] * 50.0;
      jerk_mag[k] = (varargin_1 * varargin_1 + b_varargin_1 * b_varargin_1) +
                    g_norm * g_norm;
    }
  } else {
    binary_expand_op_4(jerk_mag, r, accel_mag);
  }
  nx = jerk_mag.size(0);
  for (int k{0}; k < nx; k++) {
    jerk_mag[k] = std::sqrt(jerk_mag[k]);
  }
  features[16] = coder::mean(jerk_mag);
  features[17] = coder::b_std(jerk_mag);
  coder::b_abs(ax, r);
  coder::b_abs(ay, jerk_mag);
  coder::b_abs(az, accel_mag);
  if (r.size(0) == 1) {
    nx = jerk_mag.size(0);
  } else {
    nx = r.size(0);
  }
  if ((r.size(0) == jerk_mag.size(0)) && (nx == accel_mag.size(0))) {
    nx = r.size(0);
    for (int k{0}; k < nx; k++) {
      r[k] = (r[k] + jerk_mag[k]) + accel_mag[k];
    }
    features[18] = coder::mean(r);
  } else {
    binary_expand_op_3(features, r, jerk_mag, accel_mag);
  }
  coder::b_abs(gx, r);
  coder::b_abs(gy, jerk_mag);
  coder::b_abs(gz, accel_mag);
  if (r.size(0) == 1) {
    nx = jerk_mag.size(0);
  } else {
    nx = r.size(0);
  }
  if ((r.size(0) == jerk_mag.size(0)) && (nx == accel_mag.size(0))) {
    nx = r.size(0);
    for (int k{0}; k < nx; k++) {
      r[k] = (r[k] + jerk_mag[k]) + accel_mag[k];
    }
    features[19] = coder::mean(r);
  } else {
    binary_expand_op_2(features, r, jerk_mag, accel_mag);
  }
  //  ==================== 추가 물리 핵심 특징 ====================
  //  1) 중력 추정 및 정렬 성분
  g_vec[0] = d;
  g_vec[1] = d1;
  g_vec[2] = d2;
  g_norm = coder::b_norm(g_vec);
  if (g_norm == 0.0) {
    g_vec[0] = 0.0;
    g_vec[1] = 0.0;
    g_vec[2] = 1.0;
  } else {
    g_vec[0] = d / g_norm;
    g_vec[1] = d1 / g_norm;
    g_vec[2] = d2 / g_norm;
  }
  a.set_size(sensorData.size(0), 3);
  for (int k{0}; k < loop_ub; k++) {
    a[k] = sensorData[k];
    a[k + a.size(0)] = sensorData[k + sensorData.size(0)];
    a[k + a.size(0) * 2] = sensorData[k + sensorData.size(0) * 2];
  }
  ax.set_size(sensorData.size(0));
  if (loop_ub - 1 >= 0) {
    std::memset(&ax[0], 0, static_cast<unsigned int>(loop_ub) * sizeof(double));
  }
  for (int k{0}; k < 3; k++) {
    for (int i{0}; i < loop_ub; i++) {
      ax[i] = ax[i] + a[i + a.size(0) * k] * g_vec[k];
    }
  }
  if (a.size(0) == ax.size(0)) {
    b_a.set_size(sensorData.size(0), 3);
    for (int k{0}; k < 3; k++) {
      for (int i{0}; i < loop_ub; i++) {
        b_a[i + b_a.size(0) * k] = a[i + a.size(0) * k] - ax[i] * g_vec[k];
      }
    }
    r1.set_size(sensorData.size(0), 3);
    nx = b_a.size(0) * 3;
    for (int k{0}; k < nx; k++) {
      g_norm = b_a[k];
      r1[k] = g_norm * g_norm;
    }
    coder::sum(r1, a_perp);
  } else {
    binary_expand_op_1(a_perp, a, ax, g_vec);
  }
  nx = a_perp.size(0);
  for (int k{0}; k < nx; k++) {
    a_perp[k] = std::sqrt(a_perp[k]);
  }
  //  2) 평면성(주 운동면 비)
  g_norm = coder::var(ax);
  features[20] =
      g_norm / std::fmax(g_norm + coder::var(a_perp), 2.2204460492503131E-16);
  //  3) 손목 기울기(워치 z축 대비)
  features[21] =
      57.295779513082323 *
      std::acos(std::fmax(
          -1.0, std::fmin(1.0, (g_vec[0] * 0.0 + g_vec[1] * 0.0) + g_vec[2])));
  features[22] = coder::b_std(ax);
  //  4) 상·하 비대칭
  nx = ax.size(0);
  jerk_mag.set_size(ax.size(0));
  accel_mag.set_size(ax.size(0));
  for (int k{0}; k < nx; k++) {
    jerk_mag[k] = ax[k];
    if (ax[k] < 0.0) {
      jerk_mag[k] = 0.0;
    }
    accel_mag[k] = -ax[k];
  }
  for (int k{0}; k < nx; k++) {
    if (accel_mag[k] < 0.0) {
      accel_mag[k] = 0.0;
    }
  }
  features[23] = coder::trapz(t, jerk_mag) /
                 std::fmax(coder::trapz(t, accel_mag), 2.2204460492503131E-16);
  //  5) 전환 날카로움
  coder::gradient(ax, jerk_mag);
  nx = jerk_mag.size(0);
  for (int k{0}; k < nx; k++) {
    jerk_mag[k] = jerk_mag[k] * 50.0;
  }
  if (jerk_mag.size(0) == 0) {
    features[24] = 0.0;
  } else {
    coder::b_abs(jerk_mag, accel_mag);
    coder::internal::sort(accel_mag);
    features[24] =
        accel_mag[static_cast<int>(std::ceil(
                      0.95 * static_cast<double>(accel_mag.size(0)))) -
                  1];
  }
  //  6) 회전 휴지 비율
  nx = gyro_mag.size(0);
  b_gyro_mag.set_size(gyro_mag.size(0));
  for (int k{0}; k < nx; k++) {
    b_gyro_mag[k] = (gyro_mag[k] < 0.17453292519943295);
  }
  features[25] = coder::b_mean(b_gyro_mag);
  //  7) 이동 듀티
  coder::b_abs(ax, r);
  if (r.size(0) == gyro_mag.size(0)) {
    nx = r.size(0);
    b_gyro_mag.set_size(r.size(0));
    for (int k{0}; k < nx; k++) {
      b_gyro_mag[k] = ((r[k] > 0.15) || (gyro_mag[k] > 0.43633231299858238));
    }
    features[26] = coder::b_mean(b_gyro_mag);
  } else {
    binary_expand_op(features, r, gyro_mag);
  }
  //  8) 대략적 분당 반복수
  coder::findpeaks(ax, jerk_mag, accel_mag);
  features[27] =
      static_cast<double>(accel_mag.size(0)) /
      std::fmax(static_cast<double>(sensorData.size(0)) / 50.0 / 60.0,
                2.2204460492503131E-16);
  //  9) 경로 적분
  coder::b_abs(ax, r);
  features[28] = coder::trapz(t, r) + coder::trapz(t, a_perp);
  features[29] = coder::trapz(t, gyro_mag);
  //  10) 가속 PCA 1축성
  nx = coder::pca(a, a__4_data, a__4_size, a__5, latent_data);
  g_norm = coder::sum(latent_data, nx);
  if (g_norm == 0.0) {
    features[30] = 0.0;
  } else {
    features[30] = latent_data[0] / g_norm;
  }
  //  11) 스펙트럴 플랫니스
  coder::periodogram(ax, jerk_mag, accel_mag);
  nx = jerk_mag.size(0);
  ax.set_size(jerk_mag.size(0));
  for (int k{0}; k < nx; k++) {
    g_norm = jerk_mag[k] + 2.2204460492503131E-16;
    ax[k] = g_norm;
    jerk_mag[k] = g_norm;
  }
  features[31] = coder::geomean(ax) / coder::mean(jerk_mag);
  //  ==================== 최종 특징 벡터 할당 ====================
  //  loadLearnerForCoder로 불러온 모델은 바로 predict 함수를 사용할 수
  //  있습니다.
  modelObj.predict(features, predictedExercise);
}

void predictExercise_init()
{
  modelObj_not_empty = false;
}

// End of code generation (predictExercise.cpp)
