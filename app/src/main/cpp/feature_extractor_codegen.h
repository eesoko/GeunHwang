//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// feature_extractor_codegen.h
//
// Code generation for function 'feature_extractor_codegen'
//

#ifndef FEATURE_EXTRACTOR_CODEGEN_H
#define FEATURE_EXTRACTOR_CODEGEN_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
void binary_expand_op(double in1[32], const coder::array<double, 1U> &in2,
                      const coder::array<double, 1U> &in3);

void binary_expand_op_1(coder::array<double, 1U> &in1,
                        const coder::array<double, 2U> &in3,
                        const coder::array<double, 1U> &in4,
                        const double in5[3]);

void binary_expand_op_2(double in1[32], const coder::array<double, 1U> &in2,
                        const coder::array<double, 1U> &in3,
                        const coder::array<double, 1U> &in4);

void binary_expand_op_3(double in1[32], const coder::array<double, 1U> &in2,
                        const coder::array<double, 1U> &in3,
                        const coder::array<double, 1U> &in4);

#endif
// End of code generation (feature_extractor_codegen.h)
