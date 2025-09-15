//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// xzsvdc.h
//
// Code generation for function 'xzsvdc'
//

#ifndef XZSVDC_H
#define XZSVDC_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace internal {
namespace reflapack {
int xzsvdc(array<double, 2U> &A, array<double, 2U> &U, double S_data[],
           double V_data[], int V_size[2]);

}
} // namespace internal
} // namespace coder

#endif
// End of code generation (xzsvdc.h)
