//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// stack1.h
//
// Code generation for function 'stack1'
//

#ifndef STACK1_H
#define STACK1_H

// Include files
#include "predictExercise_internal_types.h"
#include "rtwtypes.h"
#include "coder_bounded_array.h"
#include <cstddef>
#include <cstdlib>

// Type Definitions
namespace coder {
namespace internal {
class stack {
public:
  bounded_array<struct_T, 120U, 1U> d;
  int n;
};

} // namespace internal
} // namespace coder

#endif
// End of code generation (stack1.h)
