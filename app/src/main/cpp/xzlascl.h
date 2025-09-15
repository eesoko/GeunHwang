//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// xzlascl.h
//
// Code generation for function 'xzlascl'
//

#ifndef XZLASCL_H
#define XZLASCL_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace internal {
namespace reflapack {
void xzlascl(double cfrom, double cto, int m, array<double, 2U> &A, int lda);

void xzlascl(double cfrom, double cto, int m, double A_data[]);

} // namespace reflapack
} // namespace internal
} // namespace coder

#endif
// End of code generation (xzlascl.h)
