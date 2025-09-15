//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// sum.cpp
//
// Code generation for function 'sum'
//

// Include files
#include "sum.h"
#include "blockedSummation.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
namespace coder {
void sum(const array<double, 2U> &x, array<double, 1U> &y)
{
  if (x.size(0) == 0) {
    y.set_size(0);
  } else {
    int vstride;
    vstride = x.size(0);
    y.set_size(x.size(0));
    for (int xj{0}; xj < vstride; xj++) {
      y[xj] = x[xj];
    }
    for (int xj{0}; xj < 2; xj++) {
      int xoffset;
      xoffset = (xj + 1) * vstride;
      for (int b_xj{0}; b_xj < vstride; b_xj++) {
        y[b_xj] = y[b_xj] + x[xoffset + b_xj];
      }
    }
  }
}

double sum(const double x_data[], int x_size)
{
  array<double, 1U> b_x_data;
  double y;
  if (x_size == 0) {
    y = 0.0;
  } else {
    b_x_data.set(const_cast<double *>(&x_data[0]), x_size);
    y = colMajorFlatIter(b_x_data, x_size);
  }
  return y;
}

} // namespace coder

// End of code generation (sum.cpp)
