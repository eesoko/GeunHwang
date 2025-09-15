//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// xzlangeM.cpp
//
// Code generation for function 'xzlangeM'
//

// Include files
#include "xzlangeM.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
namespace internal {
namespace reflapack {
double xzlangeM(const array<double, 2U> &x)
{
  double y;
  bool b;
  y = 0.0;
  b = (x.size(0) == 0);
  if (!b) {
    int k;
    bool exitg1;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k <= x.size(0) * 3 - 1)) {
      double absxk;
      absxk = std::abs(x[k]);
      if (std::isnan(absxk)) {
        y = rtNaN;
        exitg1 = true;
      } else {
        if (absxk > y) {
          y = absxk;
        }
        k++;
      }
    }
  }
  return y;
}

} // namespace reflapack
} // namespace internal
} // namespace coder

// End of code generation (xzlangeM.cpp)
