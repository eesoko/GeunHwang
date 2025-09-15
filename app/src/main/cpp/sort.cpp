//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// sort.cpp
//
// Code generation for function 'sort'
//

// Include files
#include "sort.h"
#include "rt_nonfinite.h"
#include "sortIdx.h"
#include "coder_array.h"
#include <cmath>
#include <cstring>

// Function Definitions
namespace coder {
namespace internal {
void sort(array<double, 1U> &x)
{
  array<double, 1U> vwork;
  array<double, 1U> xwork;
  array<int, 1U> iidx;
  array<int, 1U> iwork;
  int dim;
  int i;
  int vstride;
  dim = 2;
  if (x.size(0) != 1) {
    dim = 1;
  }
  if (dim <= 1) {
    i = x.size(0);
  } else {
    i = 1;
  }
  vwork.set_size(i);
  vstride = 1;
  dim -= 2;
  for (int k{0}; k <= dim; k++) {
    vstride *= x.size(0);
  }
  for (int j{0}; j < vstride; j++) {
    double x4[4];
    int idx4[4];
    int bLen;
    int bLen2;
    int i1;
    int i2;
    int i4;
    int ib;
    int loop_ub;
    int nBlocks;
    int wOffset;
    for (int k{0}; k < i; k++) {
      vwork[k] = x[j + k * vstride];
    }
    loop_ub = vwork.size(0);
    iidx.set_size(vwork.size(0));
    if (loop_ub - 1 >= 0) {
      std::memset(&iidx[0], 0,
                  static_cast<unsigned int>(loop_ub) * sizeof(int));
    }
    iwork.set_size(vwork.size(0));
    if (loop_ub - 1 >= 0) {
      std::memset(&iwork[0], 0,
                  static_cast<unsigned int>(loop_ub) * sizeof(int));
    }
    xwork.set_size(vwork.size(0));
    x4[0] = 0.0;
    idx4[0] = 0;
    x4[1] = 0.0;
    idx4[1] = 0;
    x4[2] = 0.0;
    idx4[2] = 0;
    x4[3] = 0.0;
    idx4[3] = 0;
    nBlocks = 0;
    ib = 0;
    for (int k{0}; k < loop_ub; k++) {
      if (std::isnan(vwork[k])) {
        dim = (loop_ub - nBlocks) - 1;
        iidx[dim] = k + 1;
        xwork[dim] = vwork[k];
        nBlocks++;
      } else {
        ib++;
        idx4[ib - 1] = k + 1;
        x4[ib - 1] = vwork[k];
        if (ib == 4) {
          double d;
          double d1;
          dim = k - nBlocks;
          if (x4[0] <= x4[1]) {
            i1 = 1;
            i2 = 2;
          } else {
            i1 = 2;
            i2 = 1;
          }
          if (x4[2] <= x4[3]) {
            ib = 3;
            i4 = 4;
          } else {
            ib = 4;
            i4 = 3;
          }
          d = x4[i1 - 1];
          d1 = x4[ib - 1];
          if (d <= d1) {
            if (x4[i2 - 1] <= d1) {
              bLen = i1;
              bLen2 = i2;
              i1 = ib;
              i2 = i4;
            } else if (x4[i2 - 1] <= x4[i4 - 1]) {
              bLen = i1;
              bLen2 = ib;
              i1 = i2;
              i2 = i4;
            } else {
              bLen = i1;
              bLen2 = ib;
              i1 = i4;
            }
          } else if (d <= x4[i4 - 1]) {
            if (x4[i2 - 1] <= x4[i4 - 1]) {
              bLen = ib;
              bLen2 = i1;
              i1 = i2;
              i2 = i4;
            } else {
              bLen = ib;
              bLen2 = i1;
              i1 = i4;
            }
          } else {
            bLen = ib;
            bLen2 = i4;
          }
          iidx[dim - 3] = idx4[bLen - 1];
          iidx[dim - 2] = idx4[bLen2 - 1];
          iidx[dim - 1] = idx4[i1 - 1];
          iidx[dim] = idx4[i2 - 1];
          vwork[dim - 3] = x4[bLen - 1];
          vwork[dim - 2] = x4[bLen2 - 1];
          vwork[dim - 1] = x4[i1 - 1];
          vwork[dim] = x4[i2 - 1];
          ib = 0;
        }
      }
    }
    wOffset = vwork.size(0) - nBlocks;
    if (ib > 0) {
      signed char perm[4];
      perm[1] = 0;
      perm[2] = 0;
      perm[3] = 0;
      if (ib == 1) {
        perm[0] = 1;
      } else if (ib == 2) {
        if (x4[0] <= x4[1]) {
          perm[0] = 1;
          perm[1] = 2;
        } else {
          perm[0] = 2;
          perm[1] = 1;
        }
      } else if (x4[0] <= x4[1]) {
        if (x4[1] <= x4[2]) {
          perm[0] = 1;
          perm[1] = 2;
          perm[2] = 3;
        } else if (x4[0] <= x4[2]) {
          perm[0] = 1;
          perm[1] = 3;
          perm[2] = 2;
        } else {
          perm[0] = 3;
          perm[1] = 1;
          perm[2] = 2;
        }
      } else if (x4[0] <= x4[2]) {
        perm[0] = 2;
        perm[1] = 1;
        perm[2] = 3;
      } else if (x4[1] <= x4[2]) {
        perm[0] = 2;
        perm[1] = 3;
        perm[2] = 1;
      } else {
        perm[0] = 3;
        perm[1] = 2;
        perm[2] = 1;
      }
      dim = static_cast<unsigned char>(ib);
      for (int k{0}; k < dim; k++) {
        i1 = (wOffset - ib) + k;
        i2 = perm[k];
        iidx[i1] = idx4[i2 - 1];
        vwork[i1] = x4[i2 - 1];
      }
    }
    dim = nBlocks >> 1;
    for (int k{0}; k < dim; k++) {
      i1 = wOffset + k;
      i2 = iidx[i1];
      ib = (loop_ub - k) - 1;
      iidx[i1] = iidx[ib];
      iidx[ib] = i2;
      vwork[i1] = xwork[ib];
      vwork[ib] = xwork[i1];
    }
    if ((static_cast<unsigned int>(nBlocks) & 1U) != 0U) {
      dim += wOffset;
      vwork[dim] = xwork[dim];
    }
    dim = 2;
    if (wOffset > 1) {
      if (vwork.size(0) >= 256) {
        nBlocks = wOffset >> 8;
        if (nBlocks > 0) {
          for (int b{0}; b < nBlocks; b++) {
            double b_xwork[256];
            int b_iwork[256];
            i4 = (b << 8) - 1;
            for (int b_b{0}; b_b < 6; b_b++) {
              bLen = 1 << (b_b + 2);
              bLen2 = bLen << 1;
              loop_ub = 256 >> (b_b + 3);
              for (int b_k{0}; b_k < loop_ub; b_k++) {
                i1 = (i4 + b_k * bLen2) + 1;
                for (int k{0}; k < bLen2; k++) {
                  dim = i1 + k;
                  b_iwork[k] = iidx[dim];
                  b_xwork[k] = vwork[dim];
                }
                i2 = 0;
                ib = bLen;
                dim = i1 - 1;
                int exitg1;
                do {
                  exitg1 = 0;
                  dim++;
                  if (b_xwork[i2] <= b_xwork[ib]) {
                    iidx[dim] = b_iwork[i2];
                    vwork[dim] = b_xwork[i2];
                    if (i2 + 1 < bLen) {
                      i2++;
                    } else {
                      exitg1 = 1;
                    }
                  } else {
                    iidx[dim] = b_iwork[ib];
                    vwork[dim] = b_xwork[ib];
                    if (ib + 1 < bLen2) {
                      ib++;
                    } else {
                      dim -= i2;
                      for (int k{i2 + 1}; k <= bLen; k++) {
                        ib = dim + k;
                        iidx[ib] = b_iwork[k - 1];
                        vwork[ib] = b_xwork[k - 1];
                      }
                      exitg1 = 1;
                    }
                  }
                } while (exitg1 == 0);
              }
            }
          }
          dim = nBlocks << 8;
          i1 = wOffset - dim;
          if (i1 > 0) {
            merge_block(iidx, vwork, dim, i1, 2, iwork, xwork);
          }
          dim = 8;
        }
      }
      merge_block(iidx, vwork, 0, wOffset, dim, iwork, xwork);
    }
    for (int k{0}; k < i; k++) {
      x[j + k * vstride] = vwork[k];
    }
  }
}

} // namespace internal
} // namespace coder

// End of code generation (sort.cpp)
