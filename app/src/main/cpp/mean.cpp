//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// mean.cpp
//
// Code generation for function 'mean'
//

// Include files
#include "mean.h"
#include "blockedSummation.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
namespace coder {
double b_mean(const array<bool, 1U> &x)
{
  double y;
  if (x.size(0) == 0) {
    y = 0.0;
  } else {
    int firstBlockLength;
    int lastBlockLength;
    int nblocks;
    if (x.size(0) <= 1024) {
      firstBlockLength = x.size(0);
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = static_cast<int>(static_cast<unsigned int>(x.size(0)) >> 10);
      lastBlockLength = x.size(0) - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    y = x[0];
    for (int k{2}; k <= firstBlockLength; k++) {
      y += static_cast<double>(x[k - 1]);
    }
    for (int k{2}; k <= nblocks; k++) {
      double bsum;
      int hi;
      firstBlockLength = (k - 1) << 10;
      bsum = x[firstBlockLength];
      if (k == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      for (int b_k{2}; b_k <= hi; b_k++) {
        bsum += static_cast<double>(x[(firstBlockLength + b_k) - 1]);
      }
      y += bsum;
    }
  }
  y /= static_cast<double>(x.size(0));
  return y;
}

double mean(const array<double, 1U> &x)
{
  double y;
  if (x.size(0) == 0) {
    y = 0.0;
  } else {
    y = colMajorFlatIter(x, x.size(0));
  }
  y /= static_cast<double>(x.size(0));
  return y;
}

} // namespace coder

// End of code generation (mean.cpp)
