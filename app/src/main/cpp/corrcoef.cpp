//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// corrcoef.cpp
//
// Code generation for function 'corrcoef'
//

// Include files
#include "corrcoef.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
void corrcoef(const array<double, 1U> &x, const array<double, 1U> &varargin_1,
              double r[4])
{
  array<double, 2U> result;
  double d[2];
  double absrij;
  double s;
  int b_i;
  int b_loop_ub;
  int loop_ub;
  loop_ub = x.size(0);
  result.set_size(x.size(0), 2);
  for (int i{0}; i < loop_ub; i++) {
    result[i] = x[i];
  }
  b_loop_ub = varargin_1.size(0);
  for (int i{0}; i < b_loop_ub; i++) {
    result[i + result.size(0)] = varargin_1[i];
  }
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  r[3] = 0.0;
  if (result.size(0) < 2) {
    r[0] = rtNaN;
    r[1] = rtNaN;
    r[2] = rtNaN;
    r[3] = rtNaN;
  } else {
    int fm;
    for (int j{0}; j < 2; j++) {
      s = 0.0;
      for (int i{0}; i < loop_ub; i++) {
        s += result[i + result.size(0) * j];
      }
      s /= static_cast<double>(result.size(0));
      for (int i{0}; i < loop_ub; i++) {
        b_loop_ub = i + result.size(0) * j;
        result[b_loop_ub] = result[b_loop_ub] - s;
      }
    }
    fm = result.size(0) - 1;
    for (int b_j{0}; b_j < 2; b_j++) {
      s = 0.0;
      for (int i{0}; i < loop_ub; i++) {
        absrij = result[i + result.size(0) * b_j];
        s += absrij * absrij;
      }
      b_loop_ub = b_j << 1;
      r[b_j + b_loop_ub] = s / static_cast<double>(fm);
      b_i = b_j + 2;
      for (int i{b_i}; i < 3; i++) {
        s = 0.0;
        for (int j{0}; j < loop_ub; j++) {
          s += result[j + result.size(0)] * result[j + result.size(0) * b_j];
        }
        r[b_loop_ub + 1] = s / static_cast<double>(fm);
      }
    }
  }
  d[0] = std::sqrt(r[0]);
  d[1] = std::sqrt(r[3]);
  for (int j{0}; j < 2; j++) {
    b_i = j + 2;
    for (int i{b_i}; i < 3; i++) {
      b_loop_ub = (j << 1) + 1;
      r[b_loop_ub] = r[b_loop_ub] / d[1] / d[j];
    }
    for (int i{b_i}; i < 3; i++) {
      b_loop_ub = (j << 1) + 1;
      s = r[b_loop_ub];
      absrij = std::abs(s);
      if (absrij > 1.0) {
        s /= absrij;
        r[b_loop_ub] = s;
      }
      r[j + 2] = r[b_loop_ub];
    }
    b_loop_ub = j + (j << 1);
    s = r[b_loop_ub];
    if (s > 0.0) {
      if (std::isnan(s)) {
        r[b_loop_ub] = rtNaN;
      } else if (s < 0.0) {
        r[b_loop_ub] = -1.0;
      } else {
        r[b_loop_ub] = (s > 0.0);
      }
    } else {
      r[b_loop_ub] = rtNaN;
    }
  }
}

} // namespace coder

// End of code generation (corrcoef.cpp)
