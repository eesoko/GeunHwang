/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * introsort.c
 *
 * Code generation for function 'introsort'
 *
 */

/* Include files */
#include "introsort.h"
#include "heapsort.h"
#include "insertionsort.h"
#include "rt_nonfinite.h"

/* Type Definitions */
#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct {
  int xstart;
  int xend;
  int depth;
} struct_T;
#endif /* typedef_struct_T */

/* Function Definitions */
void introsort(int x_data[], int xend)
{
  struct_T frame;
  int pmax;
  int pmin;
  if (xend > 1) {
    if (xend <= 32) {
      insertionsort(x_data, 1, xend);
    } else {
      struct_T st_d_data[120];
      int MAXDEPTH;
      int pivot;
      int pow2p;
      int st_n;
      boolean_T exitg1;
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
      st_d_data[0] = frame;
      st_n = 1;
      while (st_n > 0) {
        frame = st_d_data[st_n - 1];
        st_n--;
        pmax = frame.xend - frame.xstart;
        if (pmax + 1 <= 32) {
          insertionsort(x_data, frame.xstart, frame.xend);
        } else if (frame.depth == MAXDEPTH) {
          b_heapsort(x_data, frame.xstart, frame.xend);
        } else {
          int xmid;
          xmid = (frame.xstart + pmax / 2) - 1;
          pmax = x_data[frame.xstart - 1];
          if (x_data[xmid] < pmax) {
            x_data[frame.xstart - 1] = x_data[xmid];
            x_data[xmid] = pmax;
          }
          pmax = x_data[frame.xstart - 1];
          pmin = x_data[frame.xend - 1];
          if (pmin < pmax) {
            x_data[frame.xstart - 1] = pmin;
            x_data[frame.xend - 1] = pmax;
          }
          pmax = x_data[frame.xend - 1];
          if (pmax < x_data[xmid]) {
            pow2p = x_data[xmid];
            x_data[xmid] = pmax;
            x_data[frame.xend - 1] = pow2p;
          }
          pivot = x_data[xmid];
          x_data[xmid] = x_data[frame.xend - 2];
          x_data[frame.xend - 2] = pivot;
          pmax = frame.xstart - 1;
          pmin = frame.xend - 2;
          int exitg2;
          do {
            exitg2 = 0;
            for (pmax++; x_data[pmax] < pivot; pmax++) {
            }
            for (pmin--; pivot < x_data[pmin]; pmin--) {
            }
            if (pmax + 1 >= pmin + 1) {
              exitg2 = 1;
            } else {
              pow2p = x_data[pmax];
              x_data[pmax] = x_data[pmin];
              x_data[pmin] = pow2p;
            }
          } while (exitg2 == 0);
          x_data[frame.xend - 2] = x_data[pmax];
          x_data[pmax] = pivot;
          if (pmax + 2 < frame.xend) {
            st_d_data[st_n].xstart = pmax + 2;
            st_d_data[st_n].xend = frame.xend;
            st_d_data[st_n].depth = frame.depth + 1;
            st_n++;
          }
          if (frame.xstart < pmax + 1) {
            st_d_data[st_n].xstart = frame.xstart;
            st_d_data[st_n].xend = pmax + 1;
            st_d_data[st_n].depth = frame.depth + 1;
            st_n++;
          }
        }
      }
    }
  }
}

/* End of code generation (introsort.c) */
