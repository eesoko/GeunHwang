//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// introsort.cpp
//
// Code generation for function 'introsort'
//

// Include files
#include "introsort.h"
#include "anonymous_function.h"
#include "heapsort.h"
#include "insertionsort.h"
#include "predictExercise_internal_types.h"
#include "rt_nonfinite.h"
#include "stack1.h"
#include "coder_array.h"
#include "coder_bounded_array.h"

// Function Definitions
namespace coder {
namespace internal {
void introsort(array<int, 1U> &x, int xend)
{
  struct_T frame;
  int pmax;
  int pmin;
  if (xend > 1) {
    if (xend <= 32) {
      insertionsort(x, 1, xend);
    } else {
      stack st;
      int MAXDEPTH;
      int pivot;
      int pow2p;
      bool exitg1;
      pmax = 31;
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        pivot = (pmin + pmax) >> 1;
        pow2p = 1 << pivot;
        if (pow2p == xend) {
          pmax = pivot;
          exitg1 = true;
        } else if (pow2p > xend) {
          pmax = pivot;
        } else {
          pmin = pivot;
        }
      }
      MAXDEPTH = (pmax - 1) << 1;
      frame.xstart = 1;
      frame.xend = xend;
      frame.depth = 0;
      st.d.size[0] = MAXDEPTH << 1;
      st.d.data[0] = frame;
      st.n = 1;
      while (st.n > 0) {
        int frame_tmp_tmp;
        frame_tmp_tmp = st.n - 1;
        frame = st.d.data[st.n - 1];
        st.n--;
        pmax = frame.xend - frame.xstart;
        if (pmax + 1 <= 32) {
          insertionsort(x, frame.xstart, frame.xend);
        } else if (frame.depth == MAXDEPTH) {
          b_heapsort(x, frame.xstart, frame.xend);
        } else {
          int xmid;
          xmid = (frame.xstart + ((pmax + (pmax < 0)) >> 1)) - 1;
          pmax = x[frame.xstart - 1];
          if (x[xmid] < pmax) {
            x[frame.xstart - 1] = x[xmid];
            x[xmid] = pmax;
          }
          pmax = x[frame.xstart - 1];
          pmin = x[frame.xend - 1];
          if (pmin < pmax) {
            x[frame.xstart - 1] = pmin;
            x[frame.xend - 1] = pmax;
          }
          pmax = x[frame.xend - 1];
          if (pmax < x[xmid]) {
            pow2p = x[xmid];
            x[xmid] = pmax;
            x[frame.xend - 1] = pow2p;
          }
          pivot = x[xmid];
          x[xmid] = x[frame.xend - 2];
          x[frame.xend - 2] = pivot;
          pmax = frame.xstart - 1;
          pmin = frame.xend - 2;
          int exitg2;
          do {
            exitg2 = 0;
            for (pmax++; x[pmax] < pivot; pmax++) {
            }
            for (pmin--; pivot < x[pmin]; pmin--) {
            }
            if (pmax + 1 >= pmin + 1) {
              exitg2 = 1;
            } else {
              pow2p = x[pmax];
              x[pmax] = x[pmin];
              x[pmin] = pow2p;
            }
          } while (exitg2 == 0);
          x[frame.xend - 2] = x[pmax];
          x[pmax] = pivot;
          if (pmax + 2 < frame.xend) {
            st.d.data[frame_tmp_tmp].xstart = pmax + 2;
            st.d.data[frame_tmp_tmp].xend = frame.xend;
            st.d.data[frame_tmp_tmp].depth = frame.depth + 1;
            st.n = frame_tmp_tmp + 1;
          }
          if (frame.xstart < pmax + 1) {
            st.d.data[st.n].xstart = frame.xstart;
            st.d.data[st.n].xend = pmax + 1;
            st.d.data[st.n].depth = frame.depth + 1;
            st.n++;
          }
        }
      }
    }
  }
}

void introsort(const anonymous_function &cmp, int x[6])
{
  for (int i{0}; i < 6; i++) {
    x[i] = i + 1;
  }
  insertionsort(x, cmp);
}

} // namespace internal
} // namespace coder

// End of code generation (introsort.cpp)
