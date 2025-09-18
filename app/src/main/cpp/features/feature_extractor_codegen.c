/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * feature_extractor_codegen.c
 *
 * Code generation for function 'feature_extractor_codegen'
 *
 */

/* Include files */
#include "feature_extractor_codegen.h"
#include "abs.h"
#include "blockedSummation.h"
#include "corrcoef.h"
#include "feature_extractor_codegen_data.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_initialize.h"
#include "feature_extractor_codegen_types.h"
#include "findpeaks.h"
#include "geomean.h"
#include "gradient.h"
#include "mean.h"
#include "minOrMax.h"
#include "norm.h"
#include "pca.h"
#include "periodogram.h"
#include "rms.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "sqrt.h"
#include "std.h"
#include "sum.h"
#include "trapz.h"
#include "var.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Declarations */
static void binary_expand_op(double in1[32], const double in2_data[],
                             const int *in2_size, const double in3_data[],
                             const int *in3_size);

static int binary_expand_op_1(double in1_data[], const double in3_data[],
                              const int in3_size[2], const double in4_data[],
                              const int *in4_size, const double in5[3]);

static void binary_expand_op_2(double in1[32], const double in2_data[],
                               const int *in2_size, const double in3_data[],
                               const int *in3_size, const double in4_data[],
                               const int *in4_size);

static void binary_expand_op_3(double in1[32], const double in2_data[],
                               const int *in2_size, const double in3_data[],
                               const int *in3_size, const double in4_data[],
                               const int *in4_size);

static void binary_expand_op_4(double in1_data[], int *in1_size,
                               const double in3_data[], const int *in3_size,
                               double in4, const double in5_data[],
                               const int *in5_size);

/* Function Definitions */
static void binary_expand_op(double in1[32], const double in2_data[],
                             const int *in2_size, const double in3_data[],
                             const int *in3_size)
{
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  boolean_T b_in2_data[150];
  if (*in3_size == 1) {
    loop_ub = *in2_size;
  } else {
    loop_ub = *in3_size;
  }
  stride_0_0 = (*in2_size != 1);
  stride_1_0 = (*in3_size != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in2_data[i] = ((in2_data[i * stride_0_0] > 0.15) ||
                     (in3_data[i * stride_1_0] > 0.43633231299858238));
  }
  in1[26] = b_mean(b_in2_data, loop_ub);
}

static int binary_expand_op_1(double in1_data[], const double in3_data[],
                              const int in3_size[2], const double in4_data[],
                              const int *in4_size, const double in5[3])
{
  double b_in3_data[450];
  double tmp_data[450];
  int tmp_size[2];
  int i;
  int i1;
  int in1_size;
  int loop_ub;
  int stride_0_0;
  in1_size = *in4_size;
  if (in1_size == 1) {
    loop_ub = in3_size[0];
  } else {
    loop_ub = in1_size;
  }
  stride_0_0 = (in3_size[0] != 1);
  in1_size = (in1_size != 1);
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in3_data[i1 + loop_ub * i] =
          in3_data[i1 * stride_0_0 + in3_size[0] * i] -
          in4_data[i1 * in1_size] * in5[i];
    }
  }
  tmp_size[0] = loop_ub;
  tmp_size[1] = 3;
  in1_size = loop_ub * 3;
  loop_ub = in1_size / 2 * 2;
  stride_0_0 = loop_ub - 2;
  for (i = 0; i <= stride_0_0; i += 2) {
    __m128d r;
    r = _mm_loadu_pd(&b_in3_data[i]);
    r = _mm_mul_pd(r, r);
    _mm_storeu_pd(&tmp_data[i], r);
  }
  for (i = loop_ub; i < in1_size; i++) {
    double varargin_1;
    varargin_1 = b_in3_data[i];
    tmp_data[i] = varargin_1 * varargin_1;
  }
  return sum(tmp_data, tmp_size, in1_data);
}

static void binary_expand_op_2(double in1[32], const double in2_data[],
                               const int *in2_size, const double in3_data[],
                               const int *in3_size, const double in4_data[],
                               const int *in4_size)
{
  emxArray_real_T c_in2_data;
  double b_in2_data[150];
  double accumulatedData;
  int b_in2_size;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  if (*in4_size == 1) {
    if (*in3_size == 1) {
      loop_ub = *in2_size;
    } else {
      loop_ub = *in3_size;
    }
  } else {
    loop_ub = *in4_size;
  }
  b_in2_size = loop_ub;
  stride_0_0 = (*in2_size != 1);
  stride_1_0 = (*in3_size != 1);
  stride_2_0 = (*in4_size != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in2_data[i] = (in2_data[i * stride_0_0] + in3_data[i * stride_1_0]) +
                    in4_data[i * stride_2_0];
  }
  c_in2_data.data = &b_in2_data[0];
  c_in2_data.size = &b_in2_size;
  c_in2_data.allocatedSize = 150;
  c_in2_data.numDimensions = 1;
  c_in2_data.canFreeData = false;
  if ((&b_in2_size)[0] == 0) {
    accumulatedData = 0.0;
  } else {
    accumulatedData = colMajorFlatIter(&c_in2_data, (&b_in2_size)[0]);
  }
  accumulatedData /= (double)(&b_in2_size)[0];
  in1[19] = accumulatedData;
}

static void binary_expand_op_3(double in1[32], const double in2_data[],
                               const int *in2_size, const double in3_data[],
                               const int *in3_size, const double in4_data[],
                               const int *in4_size)
{
  emxArray_real_T c_in2_data;
  double b_in2_data[150];
  double accumulatedData;
  int b_in2_size;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  if (*in4_size == 1) {
    if (*in3_size == 1) {
      loop_ub = *in2_size;
    } else {
      loop_ub = *in3_size;
    }
  } else {
    loop_ub = *in4_size;
  }
  b_in2_size = loop_ub;
  stride_0_0 = (*in2_size != 1);
  stride_1_0 = (*in3_size != 1);
  stride_2_0 = (*in4_size != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in2_data[i] = (in2_data[i * stride_0_0] + in3_data[i * stride_1_0]) +
                    in4_data[i * stride_2_0];
  }
  c_in2_data.data = &b_in2_data[0];
  c_in2_data.size = &b_in2_size;
  c_in2_data.allocatedSize = 150;
  c_in2_data.numDimensions = 1;
  c_in2_data.canFreeData = false;
  if ((&b_in2_size)[0] == 0) {
    accumulatedData = 0.0;
  } else {
    accumulatedData = colMajorFlatIter(&c_in2_data, (&b_in2_size)[0]);
  }
  accumulatedData /= (double)(&b_in2_size)[0];
  in1[18] = accumulatedData;
}

static void binary_expand_op_4(double in1_data[], int *in1_size,
                               const double in3_data[], const int *in3_size,
                               double in4, const double in5_data[],
                               const int *in5_size)
{
  double tmp_data[150];
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  if (*in1_size == 1) {
    if (*in5_size == 1) {
      loop_ub = *in3_size;
    } else {
      loop_ub = *in5_size;
    }
  } else {
    loop_ub = *in1_size;
  }
  stride_0_0 = (*in3_size != 1);
  stride_1_0 = (*in5_size != 1);
  stride_2_0 = (*in1_size != 1);
  for (i = 0; i < loop_ub; i++) {
    double b_varargin_1;
    double c_varargin_1;
    double varargin_1;
    varargin_1 = in3_data[i * stride_0_0] * in4;
    b_varargin_1 = in5_data[i * stride_1_0] * in4;
    c_varargin_1 = in1_data[i * stride_2_0] * in4;
    tmp_data[i] = (varargin_1 * varargin_1 + b_varargin_1 * b_varargin_1) +
                  c_varargin_1 * c_varargin_1;
  }
  *in1_size = loop_ub;
  if (loop_ub - 1 >= 0) {
    memcpy(&in1_data[0], &tmp_data[0], (unsigned int)loop_ub * sizeof(double));
  }
}

void feature_extractor_codegen(const double raw_data_data[],
                               const int raw_data_size[2], double Fs_actual,
                               double features[32])
{
  __m128d r;
  __m128d r1;
  __m128d r2;
  emxArray_real_T b_ax_data;
  emxArray_real_T *P1;
  emxArray_real_T *f1;
  emxArray_uint8_T *y;
  double a_data[450];
  double b_a_data[450];
  double b_tmp_data[450];
  double accel_mag_data[150];
  double ax_data[150];
  double ay_data[150];
  double az_data[150];
  double gx_data[150];
  double gy_data[150];
  double gyro_mag_data[150];
  double gz_data[150];
  double jerk_mag_data[150];
  double t_data[150];
  double tmp_data[150];
  double a__4_data[9];
  double dv[4];
  double g_vec[3];
  double b_varargin_1;
  double d;
  double d1;
  double d2;
  double g_norm;
  double varargin_1;
  double *P1_data;
  double *f1_data;
  int a__4_size[2];
  int a__5_size[2];
  int a_size[2];
  int accel_mag_size;
  int ax_size;
  int ay_size;
  int az_size;
  int b_i;
  int i;
  int jerk_mag_size;
  int loop_ub;
  int t_size;
  int tmp_size;
  int vectorUB;
  unsigned char size_tmp_idx_1;
  unsigned char *y_data;
  boolean_T b_gyro_mag_data[150];
  if (!isInitialized_feature_extractor_codegen) {
    feature_extractor_codegen_initialize();
  }
  /*  =========================================================================
   */
  /*  [ SCRIPT ]    : feature_extractor_codegen.m */
  /*  [ VERSION ]   : 2.2 (Fs as Input Version) */
  /*  [ AUTHOR ]    : GeunHwang Project (Adapted by Gemini) */
  /*  [ DATE ]      : 2025-09-16 */
  /*  */
  /*  [ OVERVIEW ] */
  /*    'feature_extractor_elite.m'에서 C++ 코드 생성이 가능한 순수 계산 */
  /*    로직만을 추출한 버전. 파일 입출력, UI, 지원되지 않는 함수 제거. */
  /*    샘플링 주파수(Fs)를 동적 입력으로 받아 처리하도록 수정됨. */
  /*  =========================================================================
   */
  /*  --- 0. 입력 및 상수 정의 --- */
  loop_ub = raw_data_size[0];
  ax_size = raw_data_size[0];
  ay_size = raw_data_size[0];
  az_size = raw_data_size[0];
  if (loop_ub - 1 >= 0) {
    memcpy(&ax_data[0], &raw_data_data[0],
           (unsigned int)loop_ub * sizeof(double));
  }
  for (i = 0; i < loop_ub; i++) {
    ay_data[i] = raw_data_data[i + raw_data_size[0]];
    az_data[i] = raw_data_data[i + raw_data_size[0] * 2];
    gx_data[i] = raw_data_data[i + raw_data_size[0] * 3];
    gy_data[i] = raw_data_data[i + raw_data_size[0] * 4];
    gz_data[i] = raw_data_data[i + raw_data_size[0] * 5];
  }
  emxInit_uint8_T(&y);
  y_data = y->data;
  if (raw_data_size[0] - 1 < 0) {
    size_tmp_idx_1 = 0U;
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    size_tmp_idx_1 = (unsigned char)raw_data_size[0];
    t_size = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (unsigned char)((unsigned int)raw_data_size[0] - 1U) + 1;
    emxEnsureCapacity_uint8_T(y, t_size);
    y_data = y->data;
    t_size = (unsigned char)((unsigned int)raw_data_size[0] - 1U);
    for (i = 0; i <= t_size; i++) {
      y_data[i] = (unsigned char)i;
    }
  }
  t_size = size_tmp_idx_1;
  for (i = 0; i < t_size; i++) {
    t_data[i] = (double)y_data[i] / Fs_actual;
  }
  emxFree_uint8_T(&y);
  /*  출력 변수의 크기를 1x32로 미리 고정합니다. */
  /*  --- 기본 물리량 --- */
  accel_mag_size = raw_data_size[0];
  t_size = loop_ub / 2 * 2;
  vectorUB = t_size - 2;
  for (b_i = 0; b_i <= vectorUB; b_i += 2) {
    r = _mm_loadu_pd(&raw_data_data[b_i]);
    r1 = _mm_mul_pd(r, r);
    r = _mm_loadu_pd(&raw_data_data[b_i + raw_data_size[0]]);
    r = _mm_mul_pd(r, r);
    r1 = _mm_add_pd(r1, r);
    r = _mm_loadu_pd(&raw_data_data[b_i + raw_data_size[0] * 2]);
    r = _mm_mul_pd(r, r);
    r = _mm_add_pd(r1, r);
    _mm_storeu_pd(&accel_mag_data[b_i], r);
    r = _mm_loadu_pd(&raw_data_data[b_i + raw_data_size[0] * 3]);
    r1 = _mm_mul_pd(r, r);
    r = _mm_loadu_pd(&raw_data_data[b_i + raw_data_size[0] * 4]);
    r = _mm_mul_pd(r, r);
    r1 = _mm_add_pd(r1, r);
    r = _mm_loadu_pd(&raw_data_data[b_i + raw_data_size[0] * 5]);
    r = _mm_mul_pd(r, r);
    r = _mm_add_pd(r1, r);
    _mm_storeu_pd(&gyro_mag_data[b_i], r);
  }
  for (b_i = t_size; b_i < loop_ub; b_i++) {
    g_norm = raw_data_data[b_i];
    varargin_1 = raw_data_data[b_i + raw_data_size[0]];
    b_varargin_1 = raw_data_data[b_i + raw_data_size[0] * 2];
    accel_mag_data[b_i] = (g_norm * g_norm + varargin_1 * varargin_1) +
                          b_varargin_1 * b_varargin_1;
    g_norm = raw_data_data[b_i + raw_data_size[0] * 3];
    varargin_1 = raw_data_data[b_i + raw_data_size[0] * 4];
    b_varargin_1 = raw_data_data[b_i + raw_data_size[0] * 5];
    gyro_mag_data[b_i] = (g_norm * g_norm + varargin_1 * varargin_1) +
                         b_varargin_1 * b_varargin_1;
  }
  b_sqrt(accel_mag_data, &accel_mag_size);
  b_sqrt(gyro_mag_data, &loop_ub);
  /*  ==================== 기존 정예 특징 ==================== */
  /*  기둥 1: 자세(Orientation) */
  b_ax_data.data = &ax_data[0];
  b_ax_data.size = &ax_size;
  b_ax_data.allocatedSize = 150;
  b_ax_data.numDimensions = 1;
  b_ax_data.canFreeData = false;
  d = mean(&b_ax_data);
  features[0] = d;
  b_ax_data.data = &ay_data[0];
  b_ax_data.size = &ay_size;
  b_ax_data.allocatedSize = 150;
  b_ax_data.numDimensions = 1;
  b_ax_data.canFreeData = false;
  d1 = mean(&b_ax_data);
  features[1] = d1;
  b_ax_data.data = &az_data[0];
  b_ax_data.size = &az_size;
  b_ax_data.allocatedSize = 150;
  b_ax_data.numDimensions = 1;
  b_ax_data.canFreeData = false;
  d2 = mean(&b_ax_data);
  features[2] = d2;
  g_norm = b_std(ax_data, ax_size);
  features[3] = g_norm;
  varargin_1 = b_std(ay_data, ay_size);
  features[4] = varargin_1;
  features[5] = b_std(az_data, az_size);
  /*  기둥 2: 궤적(Trajectory) */
  features[6] = rms(gx_data, raw_data_size[0]);
  features[7] = rms(gy_data, raw_data_size[0]);
  features[8] = rms(gz_data, raw_data_size[0]);
  /*  --- 유틸 함수 (코드 생성용) --- */
  /*  Coder가 지원하는 corrcoef를 사용하고, std가 0인 경우를 방어합니다. */
  if ((g_norm == 0.0) || (varargin_1 == 0.0)) {
    features[9] = 0.0;
  } else {
    corrcoef(ax_data, ax_size, ay_data, ay_size, dv);
    features[9] = dv[2];
  }
  /*  --- 유틸 함수 (코드 생성용) --- */
  /*  Coder가 지원하는 corrcoef를 사용하고, std가 0인 경우를 방어합니다. */
  if ((varargin_1 == 0.0) || (b_std(gz_data, raw_data_size[0]) == 0.0)) {
    features[10] = 0.0;
  } else {
    corrcoef(ay_data, ay_size, gz_data, raw_data_size[0], dv);
    features[10] = dv[2];
  }
  /*  기둥 3: 주기/강도(Periodicity & Intensity) */
  features[11] = rms(accel_mag_data, accel_mag_size);
  if (accel_mag_size == 0) {
    features[12] = 0.0;
  } else {
    features[12] = maximum(accel_mag_data, accel_mag_size) -
                   minimum(accel_mag_data, accel_mag_size);
  }
  features[13] = rms(gyro_mag_data, loop_ub);
  emxInit_real_T(&P1, 1);
  emxInit_real_T(&f1, 1);
  periodogram(accel_mag_data, accel_mag_size, Fs_actual, P1, f1);
  f1_data = f1->data;
  b_maximum(P1, &vectorUB);
  features[14] = f1_data[vectorUB - 1];
  periodogram(gyro_mag_data, loop_ub, Fs_actual, P1, f1);
  f1_data = f1->data;
  b_maximum(P1, &accel_mag_size);
  features[15] = f1_data[accel_mag_size - 1];
  /*  숨겨진 보석: 운동의 질(Quality) */
  tmp_size = gradient(ax_data, ax_size, tmp_data);
  accel_mag_size = gradient(ay_data, ay_size, accel_mag_data);
  jerk_mag_size = gradient(az_data, az_size, jerk_mag_data);
  if (tmp_size == 1) {
    t_size = accel_mag_size;
  } else {
    t_size = tmp_size;
  }
  if ((tmp_size == accel_mag_size) && (t_size == jerk_mag_size)) {
    t_size = tmp_size;
    jerk_mag_size = tmp_size;
    accel_mag_size = t_size / 2 * 2;
    vectorUB = accel_mag_size - 2;
    for (i = 0; i <= vectorUB; i += 2) {
      __m128d r3;
      r = _mm_loadu_pd(&tmp_data[i]);
      r1 = _mm_mul_pd(r, _mm_set1_pd(Fs_actual));
      r = _mm_loadu_pd(&accel_mag_data[i]);
      r2 = _mm_mul_pd(r, _mm_set1_pd(Fs_actual));
      r = _mm_loadu_pd(&jerk_mag_data[i]);
      r3 = _mm_mul_pd(r, _mm_set1_pd(Fs_actual));
      r1 = _mm_mul_pd(r1, r1);
      r = _mm_mul_pd(r2, r2);
      r1 = _mm_add_pd(r1, r);
      r = _mm_mul_pd(r3, r3);
      r = _mm_add_pd(r1, r);
      _mm_storeu_pd(&jerk_mag_data[i], r);
    }
    for (i = accel_mag_size; i < t_size; i++) {
      g_norm = tmp_data[i] * Fs_actual;
      varargin_1 = accel_mag_data[i] * Fs_actual;
      b_varargin_1 = jerk_mag_data[i] * Fs_actual;
      jerk_mag_data[i] = (g_norm * g_norm + varargin_1 * varargin_1) +
                         b_varargin_1 * b_varargin_1;
    }
  } else {
    binary_expand_op_4(jerk_mag_data, &jerk_mag_size, tmp_data, &tmp_size,
                       Fs_actual, accel_mag_data, &accel_mag_size);
  }
  b_sqrt(jerk_mag_data, &jerk_mag_size);
  b_ax_data.data = &jerk_mag_data[0];
  b_ax_data.size = &jerk_mag_size;
  b_ax_data.allocatedSize = 150;
  b_ax_data.numDimensions = 1;
  b_ax_data.canFreeData = false;
  features[16] = mean(&b_ax_data);
  features[17] = b_std(jerk_mag_data, jerk_mag_size);
  tmp_size = b_abs(ax_data, ax_size, tmp_data);
  accel_mag_size = b_abs(ay_data, ay_size, accel_mag_data);
  jerk_mag_size = b_abs(az_data, az_size, jerk_mag_data);
  if (tmp_size == 1) {
    t_size = accel_mag_size;
  } else {
    t_size = tmp_size;
  }
  if ((tmp_size == accel_mag_size) && (t_size == jerk_mag_size)) {
    t_size = tmp_size;
    accel_mag_size = (tmp_size / 2) << 1;
    vectorUB = accel_mag_size - 2;
    for (i = 0; i <= vectorUB; i += 2) {
      r = _mm_loadu_pd(&tmp_data[i]);
      r1 = _mm_loadu_pd(&accel_mag_data[i]);
      r2 = _mm_loadu_pd(&jerk_mag_data[i]);
      _mm_storeu_pd(&tmp_data[i], _mm_add_pd(_mm_add_pd(r, r1), r2));
    }
    for (i = accel_mag_size; i < t_size; i++) {
      tmp_data[i] = (tmp_data[i] + accel_mag_data[i]) + jerk_mag_data[i];
    }
    b_ax_data.data = &tmp_data[0];
    b_ax_data.size = &tmp_size;
    b_ax_data.allocatedSize = 150;
    b_ax_data.numDimensions = 1;
    b_ax_data.canFreeData = false;
    features[18] = mean(&b_ax_data);
  } else {
    binary_expand_op_3(features, tmp_data, &tmp_size, accel_mag_data,
                       &accel_mag_size, jerk_mag_data, &jerk_mag_size);
  }
  tmp_size = b_abs(gx_data, raw_data_size[0], tmp_data);
  accel_mag_size = b_abs(gy_data, raw_data_size[0], accel_mag_data);
  jerk_mag_size = b_abs(gz_data, raw_data_size[0], jerk_mag_data);
  if (tmp_size == 1) {
    t_size = accel_mag_size;
  } else {
    t_size = tmp_size;
  }
  if ((tmp_size == accel_mag_size) && (t_size == jerk_mag_size)) {
    t_size = tmp_size;
    accel_mag_size = (tmp_size / 2) << 1;
    vectorUB = accel_mag_size - 2;
    for (i = 0; i <= vectorUB; i += 2) {
      r = _mm_loadu_pd(&tmp_data[i]);
      r1 = _mm_loadu_pd(&accel_mag_data[i]);
      r2 = _mm_loadu_pd(&jerk_mag_data[i]);
      _mm_storeu_pd(&tmp_data[i], _mm_add_pd(_mm_add_pd(r, r1), r2));
    }
    for (i = accel_mag_size; i < t_size; i++) {
      tmp_data[i] = (tmp_data[i] + accel_mag_data[i]) + jerk_mag_data[i];
    }
    b_ax_data.data = &tmp_data[0];
    b_ax_data.size = &tmp_size;
    b_ax_data.allocatedSize = 150;
    b_ax_data.numDimensions = 1;
    b_ax_data.canFreeData = false;
    features[19] = mean(&b_ax_data);
  } else {
    binary_expand_op_2(features, tmp_data, &tmp_size, accel_mag_data,
                       &accel_mag_size, jerk_mag_data, &jerk_mag_size);
  }
  /*  ==================== 추가 물리 핵심 특징 ==================== */
  /*  1) 중력 추정 및 정렬 성분 */
  g_vec[0] = d;
  g_vec[1] = d1;
  g_vec[2] = d2;
  g_norm = b_norm(g_vec);
  if (g_norm == 0.0) {
    g_vec[0] = 0.0;
    g_vec[1] = 0.0;
    g_vec[2] = 1.0;
  } else {
    r = _mm_loadu_pd(&g_vec[0]);
    _mm_storeu_pd(&g_vec[0], _mm_div_pd(r, _mm_set1_pd(g_norm)));
    g_vec[2] /= g_norm;
  }
  a_size[0] = raw_data_size[0];
  a_size[1] = 3;
  ax_size = raw_data_size[0];
  for (i = 0; i < loop_ub; i++) {
    a_data[i] = raw_data_data[i];
    a_data[i + a_size[0]] = raw_data_data[i + raw_data_size[0]];
    a_data[i + a_size[0] * 2] = raw_data_data[i + raw_data_size[0] * 2];
    ax_data[i] = 0.0;
  }
  accel_mag_size = (raw_data_size[0] / 2) << 1;
  t_size = accel_mag_size - 2;
  for (i = 0; i < 3; i++) {
    g_norm = g_vec[i];
    for (b_i = 0; b_i <= t_size; b_i += 2) {
      r = _mm_loadu_pd(&a_data[b_i + a_size[0] * i]);
      r1 = _mm_loadu_pd(&ax_data[b_i]);
      _mm_storeu_pd(&ax_data[b_i],
                    _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(g_norm))));
    }
    for (b_i = accel_mag_size; b_i < loop_ub; b_i++) {
      ax_data[b_i] += a_data[b_i + a_size[0] * i] * g_norm;
    }
  }
  if (raw_data_size[0] == ax_size) {
    vectorUB = accel_mag_size - 2;
    for (i = 0; i < 3; i++) {
      for (b_i = 0; b_i <= vectorUB; b_i += 2) {
        r = _mm_loadu_pd(&ax_data[b_i]);
        t_size = b_i + a_size[0] * i;
        r1 = _mm_loadu_pd(&a_data[t_size]);
        _mm_storeu_pd(&b_a_data[t_size],
                      _mm_sub_pd(r1, _mm_mul_pd(r, _mm_set1_pd(g_vec[i]))));
      }
      for (b_i = accel_mag_size; b_i < loop_ub; b_i++) {
        t_size = b_i + a_size[0] * i;
        b_a_data[t_size] = a_data[t_size] - ax_data[b_i] * g_vec[i];
      }
    }
    a_size[1] = 3;
    t_size = raw_data_size[0] * 3;
    vectorUB = t_size / 2 * 2;
    accel_mag_size = vectorUB - 2;
    for (i = 0; i <= accel_mag_size; i += 2) {
      r = _mm_loadu_pd(&b_a_data[i]);
      r = _mm_mul_pd(r, r);
      _mm_storeu_pd(&b_tmp_data[i], r);
    }
    for (i = vectorUB; i < t_size; i++) {
      g_norm = b_a_data[i];
      b_tmp_data[i] = g_norm * g_norm;
    }
    ay_size = sum(b_tmp_data, a_size, ay_data);
  } else {
    ay_size =
        binary_expand_op_1(ay_data, a_data, a_size, ax_data, &ax_size, g_vec);
  }
  b_sqrt(ay_data, &ay_size);
  /*  2) 평면성(주 운동면 비) */
  g_norm = var(ax_data, ax_size);
  features[20] =
      g_norm / fmax(g_norm + var(ay_data, ay_size), 2.2204460492503131E-16);
  /*  3) 손목 기울기(워치 z축 대비) */
  features[21] =
      57.295779513082323 *
      acos(fmax(-1.0, fmin(1.0, (g_vec[0] * 0.0 + g_vec[1] * 0.0) + g_vec[2])));
  features[22] = b_std(ax_data, ax_size);
  /*  4) 상·하 비대칭 */
  t_size = ax_size;
  accel_mag_size = ax_size;
  jerk_mag_size = ax_size;
  for (i = 0; i < t_size; i++) {
    accel_mag_data[i] = ax_data[i];
    if (ax_data[i] < 0.0) {
      accel_mag_data[i] = 0.0;
    }
    jerk_mag_data[i] = -ax_data[i];
    if (jerk_mag_data[i] < 0.0) {
      jerk_mag_data[i] = 0.0;
    }
  }
  features[23] =
      trapz(t_data, size_tmp_idx_1, accel_mag_data, accel_mag_size) /
      fmax(trapz(t_data, size_tmp_idx_1, jerk_mag_data, jerk_mag_size),
           2.2204460492503131E-16);
  /*  5) 전환 날카로움 */
  accel_mag_size = gradient(ax_data, ax_size, accel_mag_data);
  t_size = (accel_mag_size / 2) << 1;
  vectorUB = t_size - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&accel_mag_data[i]);
    _mm_storeu_pd(&accel_mag_data[i], _mm_mul_pd(r, _mm_set1_pd(Fs_actual)));
  }
  for (i = t_size; i < accel_mag_size; i++) {
    accel_mag_data[i] *= Fs_actual;
  }
  if (accel_mag_size == 0) {
    features[24] = 0.0;
  } else {
    jerk_mag_size = b_abs(accel_mag_data, accel_mag_size, jerk_mag_data);
    sort(jerk_mag_data, &jerk_mag_size);
    features[24] = jerk_mag_data[(int)ceil(0.95 * (double)jerk_mag_size) - 1];
  }
  /*  6) 회전 휴지 비율 */
  for (i = 0; i < loop_ub; i++) {
    b_gyro_mag_data[i] = (gyro_mag_data[i] < 0.17453292519943295);
  }
  features[25] = b_mean(b_gyro_mag_data, loop_ub);
  /*  7) 이동 듀티 */
  tmp_size = b_abs(ax_data, ax_size, tmp_data);
  if (tmp_size == loop_ub) {
    for (i = 0; i < tmp_size; i++) {
      b_gyro_mag_data[i] =
          ((tmp_data[i] > 0.15) || (gyro_mag_data[i] > 0.43633231299858238));
    }
    features[26] = b_mean(b_gyro_mag_data, tmp_size);
  } else {
    binary_expand_op(features, tmp_data, &tmp_size, gyro_mag_data, &loop_ub);
  }
  /*  8) 대략적 분당 반복수 */
  g_norm = 0.5 * Fs_actual;
  if (fabs(g_norm) < 4.503599627370496E+15) {
    if (g_norm >= 0.5) {
      g_norm = floor(g_norm + 0.5);
    } else if (g_norm > -0.5) {
      g_norm *= 0.0;
    } else {
      g_norm = ceil(g_norm - 0.5);
    }
  }
  findpeaks(ax_data, ax_size, g_norm, accel_mag_data, jerk_mag_data,
            &jerk_mag_size);
  features[27] =
      (double)jerk_mag_size /
      fmax((double)raw_data_size[0] / Fs_actual / 60.0, 2.2204460492503131E-16);
  /*  9) 경로 적분 */
  t_size = b_abs(ax_data, ax_size, gx_data);
  features[28] = trapz(t_data, size_tmp_idx_1, gx_data, t_size) +
                 trapz(t_data, size_tmp_idx_1, ay_data, ay_size);
  features[29] = trapz(t_data, size_tmp_idx_1, gyro_mag_data, loop_ub);
  /*  10) 가속 PCA 1축성 */
  t_size =
      pca(a_data, a_size, a__4_data, a__4_size, b_a_data, a__5_size, g_vec);
  g_norm = b_sum(g_vec, t_size);
  if (g_norm == 0.0) {
    features[30] = 0.0;
  } else {
    features[30] = g_vec[0] / g_norm;
  }
  /*  11) 스펙트럴 플랫니스 */
  periodogram(ax_data, ax_size, Fs_actual, P1, f1);
  P1_data = P1->data;
  accel_mag_size = P1->size[0];
  t_size = f1->size[0];
  f1->size[0] = P1->size[0];
  emxEnsureCapacity_real_T(f1, t_size);
  f1_data = f1->data;
  t_size = (P1->size[0] / 2) << 1;
  vectorUB = t_size - 2;
  for (b_i = 0; b_i <= vectorUB; b_i += 2) {
    r = _mm_loadu_pd(&P1_data[b_i]);
    r = _mm_add_pd(r, _mm_set1_pd(2.2204460492503131E-16));
    _mm_storeu_pd(&f1_data[b_i], r);
    _mm_storeu_pd(&P1_data[b_i], r);
  }
  for (b_i = t_size; b_i < accel_mag_size; b_i++) {
    g_norm = P1_data[b_i] + 2.2204460492503131E-16;
    f1_data[b_i] = g_norm;
    P1_data[b_i] = g_norm;
  }
  features[31] = geomean(f1) / mean(P1);
  emxFree_real_T(&f1);
  emxFree_real_T(&P1);
  /*  ==================== 최종 특징 벡터 할당 ==================== */
}

/* End of code generation (feature_extractor_codegen.c) */
