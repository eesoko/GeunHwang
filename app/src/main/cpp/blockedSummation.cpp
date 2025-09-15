//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// blockedSummation.cpp
//
// Code generation for function 'blockedSummation'
//

// Include files
#include "blockedSummation.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
namespace coder {
double colMajorFlatIter(const array<double, 1U> &x, int vlen)
{
  double y;
  int firstBlockLength;
  int lastBlockLength;
  int nblocks;
  if (vlen <= 1024) {
    firstBlockLength = vlen;
    lastBlockLength = 0;
    nblocks = 1;
  } else {
    firstBlockLength = 1024;
    nblocks = static_cast<int>(static_cast<unsigned int>(vlen) >> 10);
    lastBlockLength = vlen - (nblocks << 10);
    if (lastBlockLength > 0) {
      nblocks++;
    } else {
      lastBlockLength = 1024;
    }
  }
  y = x[0];
  for (int k{2}; k <= firstBlockLength; k++) {
    y += x[k - 1];
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
      bsum += x[(firstBlockLength + b_k) - 1];
    }
    y += bsum;
  }
  return y;
}

} // namespace coder

// End of code generation (blockedSummation.cpp)
