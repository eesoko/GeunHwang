//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// gradient.cpp
//
// Code generation for function 'gradient'
//

// Include files
#include "gradient.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
namespace coder {
void gradient(const array<double, 1U> &x, array<double, 1U> &varargout_1)
{
  int loop_ub;
  loop_ub = x.size(0);
  varargout_1.set_size(x.size(0));
  if (x.size(0) >= 2) {
    varargout_1[0] = x[1] - x[0];
    for (int j{2}; j < loop_ub; j++) {
      varargout_1[j - 1] = (x[j] - x[j - 2]) / 2.0;
    }
    varargout_1[x.size(0) - 1] = x[x.size(0) - 1] - x[x.size(0) - 2];
  } else {
    varargout_1.set_size(x.size(0));
    for (int j{0}; j < loop_ub; j++) {
      varargout_1[0] = 0.0;
    }
  }
}

} // namespace coder

// End of code generation (gradient.cpp)
