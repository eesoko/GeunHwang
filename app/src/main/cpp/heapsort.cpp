//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// heapsort.cpp
//
// Code generation for function 'heapsort'
//

// Include files
#include "heapsort.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Declarations
namespace coder {
namespace internal {
static void heapify(array<int, 1U> &x, int idx, int xstart, int xend);

}
} // namespace coder

// Function Definitions
namespace coder {
namespace internal {
static void heapify(array<int, 1U> &x, int idx, int xstart, int xend)
{
  int cmpIdx;
  int extremum;
  int extremumIdx;
  int leftIdx;
  bool changed;
  bool exitg1;
  changed = true;
  extremumIdx = (idx + xstart) - 2;
  leftIdx = ((idx << 1) + xstart) - 1;
  exitg1 = false;
  while ((!exitg1) && (leftIdx < xend)) {
    int xcmp;
    changed = false;
    extremum = x[extremumIdx];
    cmpIdx = leftIdx - 1;
    xcmp = x[leftIdx - 1];
    if (xcmp < x[leftIdx]) {
      cmpIdx = leftIdx;
      xcmp = x[leftIdx];
    }
    if (x[extremumIdx] < xcmp) {
      x[extremumIdx] = xcmp;
      x[cmpIdx] = extremum;
      extremumIdx = cmpIdx;
      leftIdx = ((((cmpIdx - xstart) + 2) << 1) + xstart) - 1;
      changed = true;
    } else {
      exitg1 = true;
    }
  }
  if (changed && (leftIdx <= xend)) {
    extremum = x[extremumIdx];
    cmpIdx = x[leftIdx - 1];
    if (x[extremumIdx] < cmpIdx) {
      x[extremumIdx] = cmpIdx;
      x[leftIdx - 1] = extremum;
    }
  }
}

void b_heapsort(array<int, 1U> &x, int xstart, int xend)
{
  int n;
  n = xend - xstart;
  for (int idx{n + 1}; idx >= 1; idx--) {
    heapify(x, idx, xstart, xend);
  }
  for (int idx{0}; idx < n; idx++) {
    int t;
    int t_tmp;
    t_tmp = (xend - idx) - 1;
    t = x[t_tmp];
    x[t_tmp] = x[xstart - 1];
    x[xstart - 1] = t;
    heapify(x, 1, xstart, t_tmp);
  }
}

} // namespace internal
} // namespace coder

// End of code generation (heapsort.cpp)
