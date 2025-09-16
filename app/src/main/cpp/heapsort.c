/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * heapsort.c
 *
 * Code generation for function 'heapsort'
 *
 */

/* Include files */
#include "heapsort.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void heapify(int x_data[], int idx, int xstart, int xend);

/* Function Definitions */
static void heapify(int x_data[], int idx, int xstart, int xend)
{
  int cmpIdx;
  int extremum;
  int extremumIdx;
  int leftIdx;
  boolean_T changed;
  boolean_T exitg1;
  changed = true;
  extremumIdx = (idx + xstart) - 2;
  leftIdx = ((idx << 1) + xstart) - 1;
  exitg1 = false;
  while ((!exitg1) && (leftIdx < xend)) {
    int xcmp;
    changed = false;
    extremum = x_data[extremumIdx];
    cmpIdx = leftIdx - 1;
    xcmp = x_data[leftIdx - 1];
    if (xcmp < x_data[leftIdx]) {
      cmpIdx = leftIdx;
      xcmp = x_data[leftIdx];
    }
    if (x_data[extremumIdx] < xcmp) {
      x_data[extremumIdx] = xcmp;
      x_data[cmpIdx] = extremum;
      extremumIdx = cmpIdx;
      leftIdx = ((((cmpIdx - xstart) + 2) << 1) + xstart) - 1;
      changed = true;
    } else {
      exitg1 = true;
    }
  }
  if (changed && (leftIdx <= xend)) {
    extremum = x_data[extremumIdx];
    cmpIdx = x_data[leftIdx - 1];
    if (x_data[extremumIdx] < cmpIdx) {
      x_data[extremumIdx] = cmpIdx;
      x_data[leftIdx - 1] = extremum;
    }
  }
}

void b_heapsort(int x_data[], int xstart, int xend)
{
  int idx;
  int n;
  n = xend - xstart;
  for (idx = n + 1; idx >= 1; idx--) {
    heapify(x_data, idx, xstart, xend);
  }
  for (idx = 0; idx < n; idx++) {
    int t;
    int t_tmp;
    t_tmp = (xend - idx) - 1;
    t = x_data[t_tmp];
    x_data[t_tmp] = x_data[xstart - 1];
    x_data[xstart - 1] = t;
    heapify(x_data, 1, xstart, t_tmp);
  }
}

/* End of code generation (heapsort.c) */
