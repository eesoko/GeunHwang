//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// xaxpy.cpp
//
// Code generation for function 'xaxpy'
//

// Include files
#include "xaxpy.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
namespace coder {
namespace internal {
namespace blas {
void xaxpy(int n, double a, const array<double, 2U> &x, int ix0,
           array<double, 1U> &y)
{
  if ((n >= 1) && (!(a == 0.0))) {
    for (int k{0}; k < n; k++) {
      y[k + 1] = y[k + 1] + a * x[(ix0 + k) - 1];
    }
  }
}

void xaxpy(int n, double a, const array<double, 1U> &x, array<double, 2U> &y,
           int iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    for (int k{0}; k < n; k++) {
      int i;
      i = (iy0 + k) - 1;
      y[i] = y[i] + a * x[k + 1];
    }
  }
}

void xaxpy(int n, double a, int ix0, array<double, 2U> &y, int iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    for (int k{0}; k < n; k++) {
      int i;
      i = (iy0 + k) - 1;
      y[i] = y[i] + a * y[(ix0 + k) - 1];
    }
  }
}

void xaxpy(double a, double y[9], int iy0)
{
  if (!(a == 0.0)) {
    y[iy0 - 1] += a * y[1];
    y[iy0] += a * y[2];
  }
}

} // namespace blas
} // namespace internal
} // namespace coder

// End of code generation (xaxpy.cpp)
