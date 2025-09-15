//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// insertionsort.h
//
// Code generation for function 'insertionsort'
//

#ifndef INSERTIONSORT_H
#define INSERTIONSORT_H

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
void insertionsort(array<int, 1U> &x, int xstart, int xend);

void insertionsort(int x[6], const anonymous_function &cmp);

} // namespace internal
} // namespace coder

#endif
// End of code generation (insertionsort.h)
