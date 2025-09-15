//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// sortIdx.h
//
// Code generation for function 'sortIdx'
//

#ifndef SORTIDX_H
#define SORTIDX_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace internal {
void merge_block(array<int, 1U> &idx, array<double, 1U> &x, int offset, int n,
                 int preSortLevel, array<int, 1U> &iwork,
                 array<double, 1U> &xwork);

}
} // namespace coder

#endif
// End of code generation (sortIdx.h)
