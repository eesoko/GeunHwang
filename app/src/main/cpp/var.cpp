//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// var.cpp
//
// Code generation for function 'var'
//

// Include files
#include "var.h"
#include "blockedSummation.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
double var(const array<double, 1U> &x)
{
  double y;
  int n;
  n = x.size(0);
  if (x.size(0) == 0) {
    y = rtNaN;
  } else if (x.size(0) == 1) {
    if ((!std::isinf(x[0])) && (!std::isnan(x[0]))) {
      y = 0.0;
    } else {
      y = rtNaN;
    }
  } else {
    double xbar;
    xbar = colMajorFlatIter(x, x.size(0)) / static_cast<double>(x.size(0));
    y = 0.0;
    for (int k{0}; k < n; k++) {
      double t;
      t = x[k] - xbar;
      y += t * t;
    }
    y /= static_cast<double>(x.size(0)) - 1.0;
  }
  return y;
}

} // namespace coder

// End of code generation (var.cpp)
