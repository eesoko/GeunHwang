//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// trapz.cpp
//
// Code generation for function 'trapz'
//

// Include files
#include "trapz.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
double trapz(const array<double, 1U> &x, const array<double, 1U> &y)
{
  array<double, 1U> c;
  double z;
  z = 0.0;
  if (y.size(0) <= 1) {
    if ((x.size(0) == 1) && (std::isinf(x[0]) || std::isnan(x[0]))) {
      z = rtNaN;
    }
  } else {
    double c1;
    int vlen;
    if (x.size(0) == 1) {
      vlen = y.size(0);
      c.set_size(y.size(0));
      for (int iac{0}; iac < vlen; iac++) {
        c[iac] = x[0];
      }
      c1 = 0.5 * x[0];
      c[0] = c1;
      c[y.size(0) - 1] = c1;
    } else {
      vlen = y.size(0);
      c.set_size(y.size(0));
      c[0] = 0.5 * (x[1] - x[0]);
      for (int iac{2}; iac < vlen; iac++) {
        c[iac - 1] = 0.5 * (x[iac] - x[iac - 2]);
      }
      c[y.size(0) - 1] = 0.5 * (x[y.size(0) - 1] - x[y.size(0) - 2]);
    }
    for (int iac{1}; vlen < 0 ? iac >= 1 : iac <= 1; iac += vlen) {
      c1 = 0.0;
      for (int ia{1}; ia <= vlen; ia++) {
        c1 += y[ia - 1] * c[ia - 1];
      }
      z += c1;
    }
  }
  return z;
}

} // namespace coder

// End of code generation (trapz.cpp)
