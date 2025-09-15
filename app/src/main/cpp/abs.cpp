//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// abs.cpp
//
// Code generation for function 'abs'
//

// Include files
#include "abs.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
void b_abs(const array<double, 1U> &x, array<double, 1U> &y)
{
  int nx;
  nx = x.size(0);
  y.set_size(x.size(0));
  for (int k{0}; k < nx; k++) {
    y[k] = std::abs(x[k]);
  }
}

} // namespace coder

// End of code generation (abs.cpp)
