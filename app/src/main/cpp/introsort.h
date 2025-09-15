//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// introsort.h
//
// Code generation for function 'introsort'
//

#ifndef INTROSORT_H
#define INTROSORT_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
namespace coder {
class anonymous_function;

}

// Function Declarations
namespace coder {
namespace internal {
void introsort(array<int, 1U> &x, int xend);

void introsort(const anonymous_function &cmp, int x[6]);

} // namespace internal
} // namespace coder

#endif
// End of code generation (introsort.h)
