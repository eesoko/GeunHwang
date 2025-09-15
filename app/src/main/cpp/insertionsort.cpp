//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// insertionsort.cpp
//
// Code generation for function 'insertionsort'
//

// Include files
#include "insertionsort.h"
#include "anonymous_function.h"
#include "predictExercise_internal_types.h"
#include "predictExercise_types.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "coder_bounded_array.h"
#include <cmath>

// Function Definitions
namespace coder {
namespace internal {
void insertionsort(array<int, 1U> &x, int xstart, int xend)
{
  int i;
  i = xstart + 1;
  for (int k{i}; k <= xend; k++) {
    int idx;
    int xc;
    bool exitg1;
    xc = x[k - 1];
    idx = k - 1;
    exitg1 = false;
    while ((!exitg1) && (idx >= xstart)) {
      int i1;
      i1 = x[idx - 1];
      if (xc < i1) {
        x[idx] = i1;
        idx--;
      } else {
        exitg1 = true;
      }
    }
    x[idx] = xc;
  }
}

void insertionsort(int x[6], const anonymous_function &cmp)
{
  for (int k{0}; k < 5; k++) {
    int idx;
    int xc;
    bool exitg1;
    xc = x[k + 1] - 1;
    idx = k;
    exitg1 = false;
    while ((!exitg1) && (idx + 1 >= 1)) {
      int b_k;
      int b_n_tmp;
      int j;
      int n;
      int n_tmp;
      bool varargout_1;
      j = x[idx];
      n_tmp = cmp.workspace.c[xc].f1.size[1];
      b_n_tmp = cmp.workspace.c[x[idx] - 1].f1.size[1];
      n = static_cast<int>(
          std::fmin(static_cast<double>(n_tmp), static_cast<double>(b_n_tmp)));
      varargout_1 = (n_tmp < b_n_tmp);
      b_k = 0;
      int exitg2;
      do {
        exitg2 = 0;
        if (b_k <= n - 1) {
          if (cmp.workspace.c[xc].f1.data[b_k] !=
              cmp.workspace.c[x[idx] - 1].f1.data[b_k]) {
            varargout_1 = (cmp.workspace.c[xc].f1.data[b_k] <
                           cmp.workspace.c[x[idx] - 1].f1.data[b_k]);
            exitg2 = 1;
          } else {
            b_k++;
          }
        } else {
          if (n_tmp == b_n_tmp) {
            varargout_1 = (xc + 1 < j);
          }
          exitg2 = 1;
        }
      } while (exitg2 == 0);
      if (varargout_1) {
        x[idx + 1] = x[idx];
        idx--;
      } else {
        exitg1 = true;
      }
    }
    x[idx + 1] = xc + 1;
  }
}

} // namespace internal
} // namespace coder

// End of code generation (insertionsort.cpp)
