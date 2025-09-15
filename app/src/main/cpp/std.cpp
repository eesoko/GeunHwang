//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// std.cpp
//
// Code generation for function 'std'
//

// Include files
#include "std.h"
#include "blockedSummation.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
double b_std(const array<double, 1U> &x)
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
    double scale;
    double xbar;
    xbar = colMajorFlatIter(x, x.size(0)) / static_cast<double>(x.size(0));
    y = 0.0;
    scale = 3.3121686421112381E-170;
    for (int k{0}; k < n; k++) {
      double d;
      d = std::abs(x[k] - xbar);
      if (d > scale) {
        double t;
        t = scale / d;
        y = y * t * t + 1.0;
        scale = d;
      } else {
        double t;
        t = d / scale;
        y += t * t;
      }
    }
    y = scale * std::sqrt(y);
    y /= std::sqrt(static_cast<double>(x.size(0)) - 1.0);
  }
  return y;
}

} // namespace coder

// End of code generation (std.cpp)
