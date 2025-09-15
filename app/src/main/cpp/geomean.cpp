//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// geomean.cpp
//
// Code generation for function 'geomean'
//

// Include files
#include "geomean.h"
#include "blockedSummation.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
double geomean(const array<double, 1U> &x)
{
  array<double, 1U> b_x;
  int loop_ub;
  loop_ub = x.size(0);
  b_x.set_size(x.size(0));
  for (int k{0}; k < loop_ub; k++) {
    b_x[k] = std::log(x[k]);
  }
  return std::exp(colMajorFlatIter(b_x, b_x.size(0)) /
                  static_cast<double>(x.size(0)));
}

} // namespace coder

// End of code generation (geomean.cpp)
