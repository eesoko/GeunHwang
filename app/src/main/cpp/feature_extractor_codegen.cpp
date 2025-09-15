//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// feature_extractor_codegen.cpp
//
// Code generation for function 'feature_extractor_codegen'
//

// Include files
#include "feature_extractor_codegen.h"
#include "blockedSummation.h"
#include "mean.h"
#include "rt_nonfinite.h"
#include "sum.h"
#include "coder_array.h"

// Function Definitions
void binary_expand_op(double in1[32], const coder::array<double, 1U> &in2,
                      const coder::array<double, 1U> &in3)
{
  coder::array<bool, 1U> b_in2;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }
  b_in2.set_size(loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in3.size(0) != 1);
  for (int i{0}; i < loop_ub; i++) {
    b_in2[i] = ((in2[i * stride_0_0] > 0.15) ||
                (in3[i * stride_1_0] > 0.43633231299858238));
  }
  in1[26] = coder::b_mean(b_in2);
}

void binary_expand_op_1(coder::array<double, 1U> &in1,
                        const coder::array<double, 2U> &in3,
                        const coder::array<double, 1U> &in4,
                        const double in5[3])
{
  coder::array<double, 2U> b_in3;
  coder::array<double, 2U> r;
  int in4_idx_0;
  int loop_ub;
  int stride_0_0;
  in4_idx_0 = in4.size(0);
  if (in4_idx_0 == 1) {
    loop_ub = in3.size(0);
  } else {
    loop_ub = in4_idx_0;
  }
  b_in3.set_size(loop_ub, 3);
  stride_0_0 = (in3.size(0) != 1);
  in4_idx_0 = (in4_idx_0 != 1);
  for (int i{0}; i < 3; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      b_in3[i1 + b_in3.size(0) * i] =
          in3[i1 * stride_0_0 + in3.size(0) * i] - in4[i1 * in4_idx_0] * in5[i];
    }
  }
  r.set_size(loop_ub, 3);
  in4_idx_0 = b_in3.size(0) * 3;
  for (int i{0}; i < in4_idx_0; i++) {
    double varargin_1;
    varargin_1 = b_in3[i];
    r[i] = varargin_1 * varargin_1;
  }
  coder::sum(r, in1);
}

void binary_expand_op_2(double in1[32], const coder::array<double, 1U> &in2,
                        const coder::array<double, 1U> &in3,
                        const coder::array<double, 1U> &in4)
{
  coder::array<double, 1U> b_in2;
  double accumulatedData;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  if (in4.size(0) == 1) {
    if (in3.size(0) == 1) {
      loop_ub = in2.size(0);
    } else {
      loop_ub = in3.size(0);
    }
  } else {
    loop_ub = in4.size(0);
  }
  b_in2.set_size(loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in3.size(0) != 1);
  stride_2_0 = (in4.size(0) != 1);
  for (int i{0}; i < loop_ub; i++) {
    b_in2[i] =
        (in2[i * stride_0_0] + in3[i * stride_1_0]) + in4[i * stride_2_0];
  }
  if (b_in2.size(0) == 0) {
    accumulatedData = 0.0;
  } else {
    accumulatedData = coder::colMajorFlatIter(b_in2, b_in2.size(0));
  }
  in1[19] = accumulatedData / static_cast<double>(b_in2.size(0));
}

void binary_expand_op_3(double in1[32], const coder::array<double, 1U> &in2,
                        const coder::array<double, 1U> &in3,
                        const coder::array<double, 1U> &in4)
{
  coder::array<double, 1U> b_in2;
  double accumulatedData;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  if (in4.size(0) == 1) {
    if (in3.size(0) == 1) {
      loop_ub = in2.size(0);
    } else {
      loop_ub = in3.size(0);
    }
  } else {
    loop_ub = in4.size(0);
  }
  b_in2.set_size(loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in3.size(0) != 1);
  stride_2_0 = (in4.size(0) != 1);
  for (int i{0}; i < loop_ub; i++) {
    b_in2[i] =
        (in2[i * stride_0_0] + in3[i * stride_1_0]) + in4[i * stride_2_0];
  }
  if (b_in2.size(0) == 0) {
    accumulatedData = 0.0;
  } else {
    accumulatedData = coder::colMajorFlatIter(b_in2, b_in2.size(0));
  }
  in1[18] = accumulatedData / static_cast<double>(b_in2.size(0));
}

// End of code generation (feature_extractor_codegen.cpp)
