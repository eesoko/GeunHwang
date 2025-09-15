//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// sortIdx.cpp
//
// Code generation for function 'sortIdx'
//

// Include files
#include "sortIdx.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Declarations
namespace coder {
namespace internal {
static void merge(array<int, 1U> &idx, array<double, 1U> &x, int offset, int np,
                  int nq, array<int, 1U> &iwork, array<double, 1U> &xwork);

}
} // namespace coder

// Function Definitions
namespace coder {
namespace internal {
static void merge(array<int, 1U> &idx, array<double, 1U> &x, int offset, int np,
                  int nq, array<int, 1U> &iwork, array<double, 1U> &xwork)
{
  if (nq != 0) {
    int iout;
    int n;
    int p;
    int q;
    n = np + nq;
    for (int j{0}; j < n; j++) {
      q = offset + j;
      iwork[j] = idx[q];
      xwork[j] = x[q];
    }
    p = 0;
    q = np;
    iout = offset - 1;
    int exitg1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork[p] <= xwork[q]) {
        idx[iout] = iwork[p];
        x[iout] = xwork[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx[iout] = iwork[q];
        x[iout] = xwork[q];
        if (q + 1 < n) {
          q++;
        } else {
          q = iout - p;
          for (int j{p + 1}; j <= np; j++) {
            iout = q + j;
            idx[iout] = iwork[j - 1];
            x[iout] = xwork[j - 1];
          }
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

void merge_block(array<int, 1U> &idx, array<double, 1U> &x, int offset, int n,
                 int preSortLevel, array<int, 1U> &iwork,
                 array<double, 1U> &xwork)
{
  int bLen;
  int nPairs;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    int tailOffset;
    if ((static_cast<unsigned int>(nPairs) & 1U) != 0U) {
      int nTail;
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }
    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (int k{0}; k < nPairs; k++) {
      merge(idx, x, offset + k * tailOffset, bLen, bLen, iwork, xwork);
    }
    bLen = tailOffset;
  }
  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

} // namespace internal
} // namespace coder

// End of code generation (sortIdx.cpp)
