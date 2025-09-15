//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// rms.cpp
//
// Code generation for function 'rms'
//

// Include files
#include "rms.h"
#include "blockedSummation.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
double rms(const array<double, 1U> &xIn)
{
  array<double, 1U> x;
  double y;
  int loop_ub;
  loop_ub = xIn.size(0);
  x.set_size(xIn.size(0));
  for (int i{0}; i < loop_ub; i++) {
    x[i] = xIn[i] * xIn[i];
  }
  if (x.size(0) == 0) {
    y = 0.0;
  } else {
    y = colMajorFlatIter(x, x.size(0));
  }
  return std::sqrt(y / static_cast<double>(x.size(0)));
}

} // namespace coder

// End of code generation (rms.cpp)
