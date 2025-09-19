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
#include "omp.h"
#include <emmintrin.h>
#include <math.h>

/* Function Declarations */
static void binary_expand_op(double in1[32], const emxArray_real_T *in2,
                             const emxArray_real_T *in3);

static void binary_expand_op_1(emxArray_real_T *in1, const emxArray_real_T *in3,
                               const emxArray_real_T *in4, const double in5[3]);

static void binary_expand_op_2(double in1[32], const emxArray_real_T *in2,
                               const emxArray_real_T *in3,
                               const emxArray_real_T *in4);

static void binary_expand_op_3(double in1[32], const emxArray_real_T *in2,
                               const emxArray_real_T *in3,
                               const emxArray_real_T *in4);

static void binary_expand_op_4(emxArray_real_T *in1, const emxArray_real_T *in3,
                               double in4, const emxArray_real_T *in5);

/* Function Definitions */
static void binary_expand_op(double in1[32], const emxArray_real_T *in2,
                             const emxArray_real_T *in3)
{
  emxArray_boolean_T *b_in2;
  const double *in2_data;
  const double *in3_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  boolean_T *b_in2_data;
  in3_data = in3->data;
  in2_data = in2->data;
  emxInit_boolean_T(&b_in2, 1);
  if (in3->size[0] == 1) {
    loop_ub = in2->size[0];
  } else {
    loop_ub = in3->size[0];
  }
  stride_0_0 = b_in2->size[0];
  b_in2->size[0] = loop_ub;
  emxEnsureCapacity_boolean_T(b_in2, stride_0_0);
  b_in2_data = b_in2->data;
  stride_0_0 = (in2->size[0] != 1);
  stride_1_0 = (in3->size[0] != 1);
  if (loop_ub < 1200) {
    for (i = 0; i < loop_ub; i++) {
      b_in2_data[i] = ((in2_data[i * stride_0_0] > 0.15) ||
                       (in3_data[i * stride_1_0] > 0.43633231299858238));
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (i = 0; i < loop_ub; i++) {
      b_in2_data[i] = ((in2_data[i * stride_0_0] > 0.15) ||
                       (in3_data[i * stride_1_0] > 0.43633231299858238));
    }
  }
  in1[26] = b_mean(b_in2);
  emxFree_boolean_T(&b_in2);
}

static void binary_expand_op_1(emxArray_real_T *in1, const emxArray_real_T *in3,
                               const emxArray_real_T *in4, const double in5[3])
{
  emxArray_real_T *b_in3;
  emxArray_real_T *r;
  const double *in3_data;
  const double *in4_data;
  double b_varargin_1;
  double *b_in3_data;
  double *r1;
  int i;
  int i1;
  int i2;
  int in4_idx_0;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  in4_data = in4->data;
  in3_data = in3->data;
  in4_idx_0 = in4->size[0];
  emxInit_real_T(&b_in3, 2);
  if (in4_idx_0 == 1) {
    loop_ub = in3->size[0];
  } else {
    loop_ub = in4_idx_0;
  }
  stride_1_0 = b_in3->size[0] * b_in3->size[1];
  b_in3->size[0] = loop_ub;
  b_in3->size[1] = 3;
  emxEnsureCapacity_real_T(b_in3, stride_1_0);
  b_in3_data = b_in3->data;
  stride_0_0 = (in3->size[0] != 1);
  stride_1_0 = (in4_idx_0 != 1);
  in4_idx_0 = (3 * loop_ub < 1200);
  if (in4_idx_0) {
    for (i = 0; i < 3; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_in3_data[i1 + b_in3->size[0] * i] =
            in3_data[i1 * stride_0_0 + in3->size[0] * i] -
            in4_data[i1 * stride_1_0] * in5[i];
      }
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(i1)

    for (i = 0; i < 3; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_in3_data[i1 + b_in3->size[0] * i] =
            in3_data[i1 * stride_0_0 + in3->size[0] * i] -
            in4_data[i1 * stride_1_0] * in5[i];
      }
    }
  }
  emxInit_real_T(&r, 2);
  stride_1_0 = r->size[0] * r->size[1];
  r->size[0] = loop_ub;
  r->size[1] = 3;
  emxEnsureCapacity_real_T(r, stride_1_0);
  r1 = r->data;
  stride_1_0 = b_in3->size[0] * 3;
  if (in4_idx_0) {
    for (i2 = 0; i2 < stride_1_0; i2++) {
      double varargin_1;
      varargin_1 = b_in3_data[i2];
      r1[i2] = varargin_1 * varargin_1;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        b_varargin_1)

    for (i2 = 0; i2 < stride_1_0; i2++) {
      b_varargin_1 = b_in3_data[i2];
      r1[i2] = b_varargin_1 * b_varargin_1;
    }
  }
  emxFree_real_T(&b_in3);
  sum(r, in1);
  emxFree_real_T(&r);
}

static void binary_expand_op_2(double in1[32], const emxArray_real_T *in2,
                               const emxArray_real_T *in3,
                               const emxArray_real_T *in4)
{
  emxArray_real_T *b_in2;
  const double *in2_data;
  const double *in3_data;
  const double *in4_data;
  double *b_in2_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  in4_data = in4->data;
  in3_data = in3->data;
  in2_data = in2->data;
  emxInit_real_T(&b_in2, 1);
  if (in4->size[0] == 1) {
    if (in3->size[0] == 1) {
      loop_ub = in2->size[0];
    } else {
      loop_ub = in3->size[0];
    }
  } else {
    loop_ub = in4->size[0];
  }
  stride_0_0 = b_in2->size[0];
  b_in2->size[0] = loop_ub;
  emxEnsureCapacity_real_T(b_in2, stride_0_0);
  b_in2_data = b_in2->data;
  stride_0_0 = (in2->size[0] != 1);
  stride_1_0 = (in3->size[0] != 1);
  stride_2_0 = (in4->size[0] != 1);
  if (loop_ub < 1200) {
    for (i = 0; i < loop_ub; i++) {
      b_in2_data[i] = (in2_data[i * stride_0_0] + in3_data[i * stride_1_0]) +
                      in4_data[i * stride_2_0];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (i = 0; i < loop_ub; i++) {
      b_in2_data[i] = (in2_data[i * stride_0_0] + in3_data[i * stride_1_0]) +
                      in4_data[i * stride_2_0];
    }
  }
  in1[19] = mean(b_in2);
  emxFree_real_T(&b_in2);
}

static void binary_expand_op_3(double in1[32], const emxArray_real_T *in2,
                               const emxArray_real_T *in3,
                               const emxArray_real_T *in4)
{
  emxArray_real_T *b_in2;
  const double *in2_data;
  const double *in3_data;
  const double *in4_data;
  double *b_in2_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  in4_data = in4->data;
  in3_data = in3->data;
  in2_data = in2->data;
  emxInit_real_T(&b_in2, 1);
  if (in4->size[0] == 1) {
    if (in3->size[0] == 1) {
      loop_ub = in2->size[0];
    } else {
      loop_ub = in3->size[0];
    }
  } else {
    loop_ub = in4->size[0];
  }
  stride_0_0 = b_in2->size[0];
  b_in2->size[0] = loop_ub;
  emxEnsureCapacity_real_T(b_in2, stride_0_0);
  b_in2_data = b_in2->data;
  stride_0_0 = (in2->size[0] != 1);
  stride_1_0 = (in3->size[0] != 1);
  stride_2_0 = (in4->size[0] != 1);
  if (loop_ub < 1200) {
    for (i = 0; i < loop_ub; i++) {
      b_in2_data[i] = (in2_data[i * stride_0_0] + in3_data[i * stride_1_0]) +
                      in4_data[i * stride_2_0];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (i = 0; i < loop_ub; i++) {
      b_in2_data[i] = (in2_data[i * stride_0_0] + in3_data[i * stride_1_0]) +
                      in4_data[i * stride_2_0];
    }
  }
  in1[18] = mean(b_in2);
  emxFree_real_T(&b_in2);
}

static void binary_expand_op_4(emxArray_real_T *in1, const emxArray_real_T *in3,
                               double in4, const emxArray_real_T *in5)
{
  emxArray_real_T *r;
  const double *in3_data;
  const double *in5_data;
  double d_varargin_1;
  double e_varargin_1;
  double f_varargin_1;
  double *in1_data;
  double *r1;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  in5_data = in5->data;
  in3_data = in3->data;
  in1_data = in1->data;
  emxInit_real_T(&r, 1);
  if (in5->size[0] == 1) {
    if (in1->size[0] == 1) {
      loop_ub = in3->size[0];
    } else {
      loop_ub = in1->size[0];
    }
  } else {
    loop_ub = in5->size[0];
  }
  stride_0_0 = r->size[0];
  r->size[0] = loop_ub;
  emxEnsureCapacity_real_T(r, stride_0_0);
  r1 = r->data;
  stride_0_0 = (in3->size[0] != 1);
  stride_1_0 = (in1->size[0] != 1);
  stride_2_0 = (in5->size[0] != 1);
  if (loop_ub < 1200) {
    for (i = 0; i < loop_ub; i++) {
      double b_varargin_1;
      double c_varargin_1;
      double varargin_1;
      varargin_1 = in3_data[i * stride_0_0] * in4;
      b_varargin_1 = in1_data[i * stride_1_0] * in4;
      c_varargin_1 = in5_data[i * stride_2_0] * in4;
      r1[i] = (varargin_1 * varargin_1 + b_varargin_1 * b_varargin_1) +
              c_varargin_1 * c_varargin_1;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        d_varargin_1, e_varargin_1, f_varargin_1)

    for (i = 0; i < loop_ub; i++) {
      d_varargin_1 = in3_data[i * stride_0_0] * in4;
      e_varargin_1 = in1_data[i * stride_1_0] * in4;
      f_varargin_1 = in5_data[i * stride_2_0] * in4;
      r1[i] = (d_varargin_1 * d_varargin_1 + e_varargin_1 * e_varargin_1) +
              f_varargin_1 * f_varargin_1;
    }
  }
  stride_0_0 = in1->size[0];
  in1->size[0] = loop_ub;
  emxEnsureCapacity_real_T(in1, stride_0_0);
  in1_data = in1->data;
  for (i1 = 0; i1 < loop_ub; i1++) {
    in1_data[i1] = r1[i1];
  }
  emxFree_real_T(&r);
}

void feature_extractor_codegen(const emxArray_real_T *raw_data,
                               double Fs_actual, double features[32])
{
  __m128d r1;
  __m128d r2;
  __m128d r3;
  __m128d r5;
  __m128d r6;
  emxArray_boolean_T *b_gyro_mag;
  emxArray_int32_T *y;
  emxArray_real_T *a;
  emxArray_real_T *a__5;
  emxArray_real_T *a_perp;
  emxArray_real_T *accel_mag;
  emxArray_real_T *ax;
  emxArray_real_T *ay;
  emxArray_real_T *az;
  emxArray_real_T *b_a;
  emxArray_real_T *gx;
  emxArray_real_T *gy;
  emxArray_real_T *gyro_mag;
  emxArray_real_T *gz;
  emxArray_real_T *jerk_mag;
  emxArray_real_T *r;
  emxArray_real_T *r4;
  emxArray_real_T *t;
  double a__4_data[9];
  double dv[4];
  double g_vec[3];
  double latent_data[3];
  const double *raw_data_data;
  double b_varargin_1;
  double c_varargin_1;
  double d;
  double d1;
  double d2;
  double d_varargin_1;
  double e_varargin_1;
  double f_varargin_1;
  double g_norm;
  double g_varargin_1;
  double h_varargin_1;
  double i_varargin_1;
  double varargin_1;
  double *ax_data;
  double *ay_data;
  double *az_data;
  double *gx_data;
  double *gy_data;
  double *gz_data;
  int a__4_size[2];
  int b_i;
  int b_loop_ub;
  int c_i;
  int c_loop_ub;
  int d_loop_ub;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int i7;
  int i8;
  int loop_ub;
  int scalarLB;
  int *y_data;
  boolean_T *gyro_mag_data;
  if (!isInitialized_feature_extractor_codegen) {
    feature_extractor_codegen_initialize();
  }
  raw_data_data = raw_data->data;
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
  /*  */
  /*    INPUT_ARGS = {coder.typeof(0, [Inf, 6], [1, 0]), coder.typeof(0)}; */
  /*    codegen feature_extractor_codegen -args INPUT_ARGS -config:lib -lang:c
   * -report */
  /*  */
  /*  =========================================================================
   */
  /*  --- 0. 입력 및 상수 정의 --- */
  emxInit_real_T(&ax, 1);
  loop_ub = raw_data->size[0];
  b_loop_ub = ax->size[0];
  ax->size[0] = raw_data->size[0];
  emxEnsureCapacity_real_T(ax, b_loop_ub);
  ax_data = ax->data;
  emxInit_real_T(&ay, 1);
  b_loop_ub = ay->size[0];
  ay->size[0] = raw_data->size[0];
  emxEnsureCapacity_real_T(ay, b_loop_ub);
  ay_data = ay->data;
  emxInit_real_T(&az, 1);
  b_loop_ub = az->size[0];
  az->size[0] = raw_data->size[0];
  emxEnsureCapacity_real_T(az, b_loop_ub);
  az_data = az->data;
  emxInit_real_T(&gx, 1);
  b_loop_ub = gx->size[0];
  gx->size[0] = raw_data->size[0];
  emxEnsureCapacity_real_T(gx, b_loop_ub);
  gx_data = gx->data;
  emxInit_real_T(&gy, 1);
  b_loop_ub = gy->size[0];
  gy->size[0] = raw_data->size[0];
  emxEnsureCapacity_real_T(gy, b_loop_ub);
  gy_data = gy->data;
  emxInit_real_T(&gz, 1);
  b_loop_ub = gz->size[0];
  gz->size[0] = raw_data->size[0];
  emxEnsureCapacity_real_T(gz, b_loop_ub);
  gz_data = gz->data;
  for (i = 0; i < loop_ub; i++) {
    ax_data[i] = raw_data_data[i];
    ay_data[i] = raw_data_data[i + raw_data->size[0]];
    az_data[i] = raw_data_data[i + raw_data->size[0] * 2];
    gx_data[i] = raw_data_data[i + raw_data->size[0] * 3];
    gy_data[i] = raw_data_data[i + raw_data->size[0] * 4];
    gz_data[i] = raw_data_data[i + raw_data->size[0] * 5];
  }
  emxInit_int32_T(&y, 2);
  y_data = y->data;
  if (raw_data->size[0] - 1 < 0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    b_loop_ub = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = raw_data->size[0];
    emxEnsureCapacity_int32_T(y, b_loop_ub);
    y_data = y->data;
    for (i = 0; i < loop_ub; i++) {
      y_data[i] = i;
    }
  }
  emxInit_real_T(&t, 1);
  c_loop_ub = y->size[1];
  b_loop_ub = t->size[0];
  t->size[0] = y->size[1];
  emxEnsureCapacity_real_T(t, b_loop_ub);
  ay_data = t->data;
  b_loop_ub = y->size[1];
  if (y->size[1] < 1200) {
    for (b_i = 0; b_i < c_loop_ub; b_i++) {
      ay_data[b_i] = (double)y_data[b_i] / Fs_actual;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      ay_data[b_i] = (double)y_data[b_i] / Fs_actual;
    }
  }
  emxFree_int32_T(&y);
  /*  출력 변수의 크기를 1x32로 미리 고정합니다. */
  /*  --- 기본 물리량 --- */
  emxInit_real_T(&accel_mag, 1);
  b_loop_ub = accel_mag->size[0];
  accel_mag->size[0] = loop_ub;
  emxEnsureCapacity_real_T(accel_mag, b_loop_ub);
  az_data = accel_mag->data;
  emxInit_real_T(&gyro_mag, 1);
  b_loop_ub = gyro_mag->size[0];
  gyro_mag->size[0] = loop_ub;
  emxEnsureCapacity_real_T(gyro_mag, b_loop_ub);
  gy_data = gyro_mag->data;
  if (raw_data->size[0] < 1200) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      g_norm = raw_data_data[i1];
      varargin_1 = raw_data_data[i1 + raw_data->size[0]];
      b_varargin_1 = raw_data_data[i1 + raw_data->size[0] * 2];
      az_data[i1] = (g_norm * g_norm + varargin_1 * varargin_1) +
                    b_varargin_1 * b_varargin_1;
      g_norm = raw_data_data[i1 + raw_data->size[0] * 3];
      varargin_1 = raw_data_data[i1 + raw_data->size[0] * 4];
      b_varargin_1 = raw_data_data[i1 + raw_data->size[0] * 5];
      gy_data[i1] = (g_norm * g_norm + varargin_1 * varargin_1) +
                    b_varargin_1 * b_varargin_1;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        c_varargin_1, d_varargin_1, e_varargin_1)

    for (i1 = 0; i1 < loop_ub; i1++) {
      c_varargin_1 = raw_data_data[i1];
      d_varargin_1 = raw_data_data[i1 + raw_data->size[0]];
      e_varargin_1 = raw_data_data[i1 + raw_data->size[0] * 2];
      az_data[i1] =
          (c_varargin_1 * c_varargin_1 + d_varargin_1 * d_varargin_1) +
          e_varargin_1 * e_varargin_1;
      c_varargin_1 = raw_data_data[i1 + raw_data->size[0] * 3];
      d_varargin_1 = raw_data_data[i1 + raw_data->size[0] * 4];
      e_varargin_1 = raw_data_data[i1 + raw_data->size[0] * 5];
      gy_data[i1] =
          (c_varargin_1 * c_varargin_1 + d_varargin_1 * d_varargin_1) +
          e_varargin_1 * e_varargin_1;
    }
  }
  b_sqrt(accel_mag);
  b_sqrt(gyro_mag);
  gy_data = gyro_mag->data;
  /*  ==================== 기존 정예 특징 ==================== */
  /*  기둥 1: 자세(Orientation) */
  d = mean(ax);
  features[0] = d;
  d1 = mean(ay);
  features[1] = d1;
  d2 = mean(az);
  features[2] = d2;
  g_norm = b_std(ax);
  features[3] = g_norm;
  varargin_1 = b_std(ay);
  features[4] = varargin_1;
  features[5] = b_std(az);
  /*  기둥 2: 궤적(Trajectory) */
  features[6] = rms(gx);
  features[7] = rms(gy);
  features[8] = rms(gz);
  /*  --- 유틸 함수 (코드 생성용) --- */
  /*  Coder가 지원하는 corrcoef를 사용하고, std가 0인 경우를 방어합니다. */
  if ((g_norm == 0.0) || (varargin_1 == 0.0)) {
    features[9] = 0.0;
  } else {
    corrcoef(ax, ay, dv);
    features[9] = dv[2];
  }
  /*  --- 유틸 함수 (코드 생성용) --- */
  /*  Coder가 지원하는 corrcoef를 사용하고, std가 0인 경우를 방어합니다. */
  if ((varargin_1 == 0.0) || (b_std(gz) == 0.0)) {
    features[10] = 0.0;
  } else {
    corrcoef(ay, gz, dv);
    features[10] = dv[2];
  }
  /*  기둥 3: 주기/강도(Periodicity & Intensity) */
  features[11] = rms(accel_mag);
  if (accel_mag->size[0] == 0) {
    features[12] = 0.0;
  } else {
    features[12] = maximum(accel_mag) - minimum(accel_mag);
  }
  features[13] = rms(gyro_mag);
  emxInit_real_T(&jerk_mag, 1);
  emxInit_real_T(&a_perp, 1);
  periodogram(accel_mag, Fs_actual, jerk_mag, a_perp);
  ay_data = a_perp->data;
  b_maximum(jerk_mag, &b_loop_ub);
  features[14] = ay_data[b_loop_ub - 1];
  periodogram(gyro_mag, Fs_actual, accel_mag, jerk_mag);
  gx_data = jerk_mag->data;
  b_maximum(accel_mag, &b_loop_ub);
  features[15] = gx_data[b_loop_ub - 1];
  /*  숨겨진 보석: 운동의 질(Quality) */
  emxInit_real_T(&r, 1);
  gradient(ax, r);
  ay_data = r->data;
  gradient(ay, jerk_mag);
  gradient(az, accel_mag);
  az_data = accel_mag->data;
  if (r->size[0] == 1) {
    b_loop_ub = jerk_mag->size[0];
  } else {
    b_loop_ub = r->size[0];
  }
  if ((r->size[0] == jerk_mag->size[0]) && (b_loop_ub == accel_mag->size[0])) {
    d_loop_ub = r->size[0];
    b_loop_ub = jerk_mag->size[0];
    jerk_mag->size[0] = r->size[0];
    emxEnsureCapacity_real_T(jerk_mag, b_loop_ub);
    gx_data = jerk_mag->data;
    b_loop_ub = r->size[0];
    if (r->size[0] < 1200) {
      for (i2 = 0; i2 < d_loop_ub; i2++) {
        g_norm = ay_data[i2] * Fs_actual;
        varargin_1 = gx_data[i2] * Fs_actual;
        b_varargin_1 = az_data[i2] * Fs_actual;
        gx_data[i2] = (g_norm * g_norm + varargin_1 * varargin_1) +
                      b_varargin_1 * b_varargin_1;
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        f_varargin_1, g_varargin_1, h_varargin_1)

      for (i2 = 0; i2 < b_loop_ub; i2++) {
        f_varargin_1 = ay_data[i2] * Fs_actual;
        g_varargin_1 = gx_data[i2] * Fs_actual;
        h_varargin_1 = az_data[i2] * Fs_actual;
        gx_data[i2] =
            (f_varargin_1 * f_varargin_1 + g_varargin_1 * g_varargin_1) +
            h_varargin_1 * h_varargin_1;
      }
    }
  } else {
    binary_expand_op_4(jerk_mag, r, Fs_actual, accel_mag);
  }
  b_sqrt(jerk_mag);
  features[16] = mean(jerk_mag);
  features[17] = b_std(jerk_mag);
  b_abs(ax, r);
  ay_data = r->data;
  b_abs(ay, jerk_mag);
  gx_data = jerk_mag->data;
  emxFree_real_T(&ay);
  b_abs(az, accel_mag);
  az_data = accel_mag->data;
  emxFree_real_T(&az);
  if (r->size[0] == 1) {
    b_loop_ub = jerk_mag->size[0];
  } else {
    b_loop_ub = r->size[0];
  }
  if ((r->size[0] == jerk_mag->size[0]) && (b_loop_ub == accel_mag->size[0])) {
    b_loop_ub = r->size[0];
    c_loop_ub = (r->size[0] / 2) << 1;
    d_loop_ub = c_loop_ub - 2;
    for (i = 0; i <= d_loop_ub; i += 2) {
      r1 = _mm_loadu_pd(&ay_data[i]);
      r2 = _mm_loadu_pd(&gx_data[i]);
      r3 = _mm_loadu_pd(&az_data[i]);
      _mm_storeu_pd(&ay_data[i], _mm_add_pd(_mm_add_pd(r1, r2), r3));
    }
    for (i = c_loop_ub; i < b_loop_ub; i++) {
      ay_data[i] = (ay_data[i] + gx_data[i]) + az_data[i];
    }
    features[18] = mean(r);
  } else {
    binary_expand_op_3(features, r, jerk_mag, accel_mag);
  }
  b_abs(gx, r);
  ay_data = r->data;
  emxFree_real_T(&gx);
  b_abs(gy, jerk_mag);
  gx_data = jerk_mag->data;
  emxFree_real_T(&gy);
  b_abs(gz, accel_mag);
  az_data = accel_mag->data;
  emxFree_real_T(&gz);
  if (r->size[0] == 1) {
    b_loop_ub = jerk_mag->size[0];
  } else {
    b_loop_ub = r->size[0];
  }
  if ((r->size[0] == jerk_mag->size[0]) && (b_loop_ub == accel_mag->size[0])) {
    b_loop_ub = r->size[0];
    c_loop_ub = (r->size[0] / 2) << 1;
    d_loop_ub = c_loop_ub - 2;
    for (i = 0; i <= d_loop_ub; i += 2) {
      r1 = _mm_loadu_pd(&ay_data[i]);
      r2 = _mm_loadu_pd(&gx_data[i]);
      r3 = _mm_loadu_pd(&az_data[i]);
      _mm_storeu_pd(&ay_data[i], _mm_add_pd(_mm_add_pd(r1, r2), r3));
    }
    for (i = c_loop_ub; i < b_loop_ub; i++) {
      ay_data[i] = (ay_data[i] + gx_data[i]) + az_data[i];
    }
    features[19] = mean(r);
  } else {
    binary_expand_op_2(features, r, jerk_mag, accel_mag);
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
    r1 = _mm_loadu_pd(&g_vec[0]);
    _mm_storeu_pd(&g_vec[0], _mm_div_pd(r1, _mm_set1_pd(g_norm)));
    g_vec[2] /= g_norm;
  }
  emxInit_real_T(&a, 2);
  b_loop_ub = a->size[0] * a->size[1];
  a->size[0] = loop_ub;
  a->size[1] = 3;
  emxEnsureCapacity_real_T(a, b_loop_ub);
  ay_data = a->data;
  for (i = 0; i < loop_ub; i++) {
    ay_data[i] = raw_data_data[i];
    ay_data[i + a->size[0]] = raw_data_data[i + raw_data->size[0]];
    ay_data[i + a->size[0] * 2] = raw_data_data[i + raw_data->size[0] * 2];
  }
  b_loop_ub = ax->size[0];
  ax->size[0] = loop_ub;
  emxEnsureCapacity_real_T(ax, b_loop_ub);
  ax_data = ax->data;
  for (i = 0; i < loop_ub; i++) {
    ax_data[i] = 0.0;
  }
  scalarLB = (a->size[0] / 2) << 1;
  b_loop_ub = scalarLB - 2;
  for (i = 0; i < 3; i++) {
    g_norm = g_vec[i];
    for (i3 = 0; i3 <= b_loop_ub; i3 += 2) {
      r1 = _mm_loadu_pd(&ay_data[i3 + a->size[0] * i]);
      r2 = _mm_loadu_pd(&ax_data[i3]);
      _mm_storeu_pd(&ax_data[i3],
                    _mm_add_pd(r2, _mm_mul_pd(r1, _mm_set1_pd(g_norm))));
    }
    for (i3 = scalarLB; i3 < loop_ub; i3++) {
      ax_data[i3] += ay_data[i3 + a->size[0] * i] * g_norm;
    }
  }
  if (a->size[0] == ax->size[0]) {
    emxInit_real_T(&b_a, 2);
    b_loop_ub = b_a->size[0] * b_a->size[1];
    b_a->size[0] = loop_ub;
    b_a->size[1] = 3;
    emxEnsureCapacity_real_T(b_a, b_loop_ub);
    az_data = b_a->data;
    d_loop_ub = scalarLB - 2;
    b_loop_ub = (scalarLB - 1) / 2;
    c_loop_ub = a->size[0] - scalarLB;
    if (b_loop_ub >= c_loop_ub) {
      c_loop_ub = b_loop_ub;
    }
    if (3 * c_loop_ub < 1200) {
      for (i4 = 0; i4 < 3; i4++) {
        for (i5 = 0; i5 <= d_loop_ub; i5 += 2) {
          r5 = _mm_loadu_pd(&ax_data[i5]);
          r6 = _mm_loadu_pd(&ay_data[i5 + a->size[0] * i4]);
          _mm_storeu_pd(&az_data[i5 + b_a->size[0] * i4],
                        _mm_sub_pd(r6, _mm_mul_pd(r5, _mm_set1_pd(g_vec[i4]))));
        }
        for (i5 = scalarLB; i5 < loop_ub; i5++) {
          az_data[i5 + b_a->size[0] * i4] =
              ay_data[i5 + a->size[0] * i4] - ax_data[i5] * g_vec[i4];
        }
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(r5, r6, i5)

      for (i4 = 0; i4 < 3; i4++) {
        for (i5 = 0; i5 <= d_loop_ub; i5 += 2) {
          r5 = _mm_loadu_pd(&ax_data[i5]);
          r5 = _mm_mul_pd(r5, _mm_set1_pd(g_vec[i4]));
          r6 = _mm_loadu_pd(&ay_data[i5 + a->size[0] * i4]);
          r5 = _mm_sub_pd(r6, r5);
          _mm_storeu_pd(&az_data[i5 + b_a->size[0] * i4], r5);
        }
        for (i5 = scalarLB; i5 < loop_ub; i5++) {
          az_data[i5 + b_a->size[0] * i4] =
              ay_data[i5 + a->size[0] * i4] - ax_data[i5] * g_vec[i4];
        }
      }
    }
    emxInit_real_T(&r4, 2);
    b_loop_ub = r4->size[0] * r4->size[1];
    r4->size[0] = b_a->size[0];
    r4->size[1] = 3;
    emxEnsureCapacity_real_T(r4, b_loop_ub);
    ay_data = r4->data;
    b_loop_ub = b_a->size[0] * 3;
    if (b_loop_ub < 1200) {
      for (i6 = 0; i6 < b_loop_ub; i6++) {
        g_norm = az_data[i6];
        ay_data[i6] = g_norm * g_norm;
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        i_varargin_1)

      for (i6 = 0; i6 < b_loop_ub; i6++) {
        i_varargin_1 = az_data[i6];
        ay_data[i6] = i_varargin_1 * i_varargin_1;
      }
    }
    emxFree_real_T(&b_a);
    sum(r4, a_perp);
    emxFree_real_T(&r4);
  } else {
    binary_expand_op_1(a_perp, a, ax, g_vec);
  }
  b_sqrt(a_perp);
  /*  2) 평면성(주 운동면 비) */
  g_norm = var(ax);
  features[20] = g_norm / fmax(g_norm + var(a_perp), 2.2204460492503131E-16);
  /*  3) 손목 기울기(워치 z축 대비) */
  features[21] =
      57.295779513082323 *
      acos(fmax(-1.0, fmin(1.0, (g_vec[0] * 0.0 + g_vec[1] * 0.0) + g_vec[2])));
  features[22] = b_std(ax);
  /*  4) 상·하 비대칭 */
  d_loop_ub = ax->size[0];
  b_loop_ub = jerk_mag->size[0];
  jerk_mag->size[0] = ax->size[0];
  emxEnsureCapacity_real_T(jerk_mag, b_loop_ub);
  gx_data = jerk_mag->data;
  c_loop_ub = ax->size[0];
  b_loop_ub = accel_mag->size[0];
  accel_mag->size[0] = ax->size[0];
  emxEnsureCapacity_real_T(accel_mag, b_loop_ub);
  az_data = accel_mag->data;
  if (ax->size[0] < 1200) {
    for (c_i = 0; c_i < d_loop_ub; c_i++) {
      gx_data[c_i] = ax_data[c_i];
      if (ax_data[c_i] < 0.0) {
        gx_data[c_i] = 0.0;
      }
      az_data[c_i] = -ax_data[c_i];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (c_i = 0; c_i < c_loop_ub; c_i++) {
      gx_data[c_i] = ax_data[c_i];
      if (ax_data[c_i] < 0.0) {
        gx_data[c_i] = 0.0;
      }
      az_data[c_i] = -ax_data[c_i];
    }
  }
  for (i = 0; i < d_loop_ub; i++) {
    if (az_data[i] < 0.0) {
      az_data[i] = 0.0;
    }
  }
  features[23] =
      trapz(t, jerk_mag) / fmax(trapz(t, accel_mag), 2.2204460492503131E-16);
  /*  5) 전환 날카로움 */
  gradient(ax, jerk_mag);
  gx_data = jerk_mag->data;
  b_loop_ub = jerk_mag->size[0];
  c_loop_ub = (jerk_mag->size[0] / 2) << 1;
  d_loop_ub = c_loop_ub - 2;
  for (i = 0; i <= d_loop_ub; i += 2) {
    r1 = _mm_loadu_pd(&gx_data[i]);
    _mm_storeu_pd(&gx_data[i], _mm_mul_pd(r1, _mm_set1_pd(Fs_actual)));
  }
  for (i = c_loop_ub; i < b_loop_ub; i++) {
    gx_data[i] *= Fs_actual;
  }
  if (jerk_mag->size[0] == 0) {
    features[24] = 0.0;
  } else {
    b_abs(jerk_mag, accel_mag);
    sort(accel_mag);
    az_data = accel_mag->data;
    features[24] = az_data[(int)ceil(0.95 * (double)accel_mag->size[0]) - 1];
  }
  /*  6) 회전 휴지 비율 */
  emxInit_boolean_T(&b_gyro_mag, 1);
  c_loop_ub = gyro_mag->size[0];
  b_loop_ub = b_gyro_mag->size[0];
  b_gyro_mag->size[0] = gyro_mag->size[0];
  emxEnsureCapacity_boolean_T(b_gyro_mag, b_loop_ub);
  gyro_mag_data = b_gyro_mag->data;
  b_loop_ub = gyro_mag->size[0];
  if (gyro_mag->size[0] < 1200) {
    for (i7 = 0; i7 < c_loop_ub; i7++) {
      gyro_mag_data[i7] = (gy_data[i7] < 0.17453292519943295);
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (i7 = 0; i7 < b_loop_ub; i7++) {
      gyro_mag_data[i7] = (gy_data[i7] < 0.17453292519943295);
    }
  }
  features[25] = b_mean(b_gyro_mag);
  /*  7) 이동 듀티 */
  b_abs(ax, r);
  ay_data = r->data;
  if (r->size[0] == gyro_mag->size[0]) {
    c_loop_ub = r->size[0];
    b_loop_ub = b_gyro_mag->size[0];
    b_gyro_mag->size[0] = r->size[0];
    emxEnsureCapacity_boolean_T(b_gyro_mag, b_loop_ub);
    gyro_mag_data = b_gyro_mag->data;
    b_loop_ub = r->size[0];
    if (r->size[0] < 1200) {
      for (i8 = 0; i8 < c_loop_ub; i8++) {
        gyro_mag_data[i8] =
            ((ay_data[i8] > 0.15) || (gy_data[i8] > 0.43633231299858238));
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

      for (i8 = 0; i8 < b_loop_ub; i8++) {
        gyro_mag_data[i8] =
            ((ay_data[i8] > 0.15) || (gy_data[i8] > 0.43633231299858238));
      }
    }
    features[26] = b_mean(b_gyro_mag);
  } else {
    binary_expand_op(features, r, gyro_mag);
  }
  emxFree_boolean_T(&b_gyro_mag);
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
  findpeaks(ax, g_norm, jerk_mag, accel_mag);
  features[27] = (double)accel_mag->size[0] /
                 fmax((double)raw_data->size[0] / Fs_actual / 60.0,
                      2.2204460492503131E-16);
  /*  9) 경로 적분 */
  b_abs(ax, r);
  features[28] = trapz(t, r) + trapz(t, a_perp);
  emxFree_real_T(&r);
  emxFree_real_T(&a_perp);
  features[29] = trapz(t, gyro_mag);
  emxFree_real_T(&gyro_mag);
  emxFree_real_T(&t);
  /*  10) 가속 PCA 1축성 */
  emxInit_real_T(&a__5, 2);
  b_loop_ub = pca(a, a__4_data, a__4_size, a__5, latent_data);
  emxFree_real_T(&a__5);
  emxFree_real_T(&a);
  g_norm = b_sum(latent_data, b_loop_ub);
  if (g_norm == 0.0) {
    features[30] = 0.0;
  } else {
    features[30] = latent_data[0] / g_norm;
  }
  /*  11) 스펙트럴 플랫니스 */
  periodogram(ax, Fs_actual, jerk_mag, accel_mag);
  gx_data = jerk_mag->data;
  emxFree_real_T(&accel_mag);
  d_loop_ub = jerk_mag->size[0];
  b_loop_ub = ax->size[0];
  ax->size[0] = jerk_mag->size[0];
  emxEnsureCapacity_real_T(ax, b_loop_ub);
  ax_data = ax->data;
  b_loop_ub = (jerk_mag->size[0] / 2) << 1;
  c_loop_ub = b_loop_ub - 2;
  for (i = 0; i <= c_loop_ub; i += 2) {
    r1 = _mm_loadu_pd(&gx_data[i]);
    r1 = _mm_add_pd(r1, _mm_set1_pd(2.2204460492503131E-16));
    _mm_storeu_pd(&ax_data[i], r1);
    _mm_storeu_pd(&gx_data[i], r1);
  }
  for (i = b_loop_ub; i < d_loop_ub; i++) {
    g_norm = gx_data[i] + 2.2204460492503131E-16;
    ax_data[i] = g_norm;
    gx_data[i] = g_norm;
  }
  features[31] = geomean(ax) / mean(jerk_mag);
  emxFree_real_T(&jerk_mag);
  emxFree_real_T(&ax);
  /*  ==================== 최종 특징 벡터 할당 ==================== */
}

/* End of code generation (feature_extractor_codegen.c) */
