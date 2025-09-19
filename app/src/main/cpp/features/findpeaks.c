/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * findpeaks.c
 *
 * Code generation for function 'findpeaks'
 *
 */

/* Include files */
#include "findpeaks.h"
#include "eml_setop.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "introsort.h"
#include "rt_nonfinite.h"
#include "omp.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void findpeaks(const emxArray_real_T *Yin, double varargin_2,
               emxArray_real_T *Ypk, emxArray_real_T *Xpk)
{
  emxArray_boolean_T *idelete;
  emxArray_int32_T *c;
  emxArray_int32_T *iInfinite;
  emxArray_int32_T *iPk;
  emxArray_int32_T *idx;
  emxArray_int32_T *iwork;
  emxArray_int32_T *r;
  emxArray_int32_T *sortIdx;
  emxArray_int32_T *y;
  emxArray_uint32_T *x;
  const double *Yin_data;
  double ykfirst;
  double *Xpk_data;
  double *Ypk_data;
  int b_k;
  int i;
  int k;
  int kfirst;
  int loop_ub;
  int nInf;
  int nPk;
  int qEnd;
  int *c_data;
  int *iPk_data;
  int *idx_data;
  int *sortIdx_data;
  unsigned int *x_data;
  char dir;
  boolean_T isinfykfirst;
  boolean_T *idelete_data;
  Yin_data = Yin->data;
  emxInit_uint32_T(&x);
  loop_ub = Yin->size[0];
  kfirst = x->size[0];
  x->size[0] = Yin->size[0];
  emxEnsureCapacity_uint32_T(x, kfirst);
  x_data = x->data;
  kfirst = Yin->size[0];
  if (Yin->size[0] < 1200) {
    for (i = 0; i < loop_ub; i++) {
      x_data[i] = (unsigned int)i + 1U;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (i = 0; i < kfirst; i++) {
      x_data[i] = (unsigned int)i + 1U;
    }
  }
  emxInit_int32_T(&idx, 1);
  kfirst = idx->size[0];
  idx->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(idx, kfirst);
  idx_data = idx->data;
  emxInit_int32_T(&iInfinite, 1);
  kfirst = iInfinite->size[0];
  iInfinite->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(iInfinite, kfirst);
  iPk_data = iInfinite->data;
  emxInit_int32_T(&sortIdx, 1);
  nPk = 0;
  nInf = 0;
  dir = 'n';
  kfirst = 0;
  ykfirst = rtInf;
  isinfykfirst = true;
  for (k = 1; k <= loop_ub; k++) {
    double yk;
    boolean_T isinfyk;
    yk = Yin_data[k - 1];
    if (rtIsNaN(yk)) {
      yk = rtInf;
      isinfyk = true;
    } else if (rtIsInf(yk) && (yk > 0.0)) {
      isinfyk = true;
      nInf++;
      iPk_data[nInf - 1] = k;
    } else {
      isinfyk = false;
    }
    if (yk != ykfirst) {
      char previousdir;
      previousdir = dir;
      if (isinfyk || isinfykfirst) {
        dir = 'n';
      } else if (yk < ykfirst) {
        dir = 'd';
        if (previousdir == 'i') {
          nPk++;
          idx_data[nPk - 1] = kfirst;
        }
      } else {
        dir = 'i';
      }
      ykfirst = yk;
      kfirst = k;
      isinfykfirst = isinfyk;
    }
  }
  if (nPk < 1) {
    nPk = 0;
  }
  kfirst = idx->size[0];
  idx->size[0] = nPk;
  emxEnsureCapacity_int32_T(idx, kfirst);
  idx_data = idx->data;
  kfirst = iInfinite->size[0];
  if (nInf < 1) {
    iInfinite->size[0] = 0;
  } else {
    iInfinite->size[0] = nInf;
  }
  emxEnsureCapacity_int32_T(iInfinite, kfirst);
  emxInit_int32_T(&iPk, 1);
  kfirst = iPk->size[0];
  iPk->size[0] = nPk;
  emxEnsureCapacity_int32_T(iPk, kfirst);
  iPk_data = iPk->data;
  nInf = 0;
  for (k = 0; k < nPk; k++) {
    ykfirst = Yin_data[idx_data[k] - 1];
    if ((ykfirst > rtMinusInf) &&
        (ykfirst - fmax(Yin_data[idx_data[k] - 2], Yin_data[idx_data[k]]) >=
         0.0)) {
      nInf++;
      iPk_data[nInf - 1] = idx_data[k];
    }
  }
  kfirst = iPk->size[0];
  iPk->size[0] = nInf;
  emxEnsureCapacity_int32_T(iPk, kfirst);
  emxInit_int32_T(&c, 1);
  do_vectors(iPk, iInfinite, c, sortIdx, idx);
  c_data = c->data;
  emxFree_int32_T(&iInfinite);
  emxFree_int32_T(&iPk);
  emxInit_boolean_T(&idelete, 1);
  emxInit_int32_T(&r, 1);
  emxInit_int32_T(&y, 2);
  emxInit_int32_T(&iwork, 1);
  if ((c->size[0] == 0) || (varargin_2 == 0.0)) {
    nInf = c->size[0];
    kfirst = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = c->size[0];
    emxEnsureCapacity_int32_T(y, kfirst);
    iPk_data = y->data;
    if (c->size[0] > 0) {
      iPk_data[0] = 1;
      kfirst = 1;
      for (k = 2; k <= nInf; k++) {
        kfirst++;
        iPk_data[k - 1] = kfirst;
      }
    }
    kfirst = idx->size[0];
    idx->size[0] = c->size[0];
    emxEnsureCapacity_int32_T(idx, kfirst);
    idx_data = idx->data;
    for (k = 0; k < nInf; k++) {
      idx_data[k] = iPk_data[k];
    }
  } else {
    int b_i;
    int b_loop_ub;
    int n;
    b_loop_ub = c->size[0];
    kfirst = sortIdx->size[0];
    sortIdx->size[0] = c->size[0];
    emxEnsureCapacity_int32_T(sortIdx, kfirst);
    sortIdx_data = sortIdx->data;
    n = c->size[0] + 1;
    kfirst = iwork->size[0];
    iwork->size[0] = c->size[0];
    emxEnsureCapacity_int32_T(iwork, kfirst);
    iPk_data = iwork->data;
    kfirst = c->size[0] - 1;
    for (k = 1; k <= kfirst; k += 2) {
      ykfirst = Yin_data[c_data[k - 1] - 1];
      if ((ykfirst >= Yin_data[c_data[k] - 1]) || rtIsNaN(ykfirst)) {
        sortIdx_data[k - 1] = k;
        sortIdx_data[k] = k + 1;
      } else {
        sortIdx_data[k - 1] = k + 1;
        sortIdx_data[k] = k;
      }
    }
    if (((unsigned int)c->size[0] & 1U) != 0U) {
      sortIdx_data[c->size[0] - 1] = c->size[0];
    }
    b_i = 2;
    while (b_i < n - 1) {
      int i2;
      int j;
      i2 = b_i << 1;
      j = 1;
      for (kfirst = b_i + 1; kfirst < n; kfirst = qEnd + b_i) {
        int c_k;
        int kEnd;
        nInf = j - 1;
        nPk = kfirst;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }
        c_k = 0;
        kEnd = qEnd - j;
        while (c_k < kEnd) {
          int i1;
          ykfirst = Yin_data[c_data[sortIdx_data[nInf] - 1] - 1];
          i1 = sortIdx_data[nPk - 1];
          if ((ykfirst >= Yin_data[c_data[i1 - 1] - 1]) || rtIsNaN(ykfirst)) {
            iPk_data[c_k] = sortIdx_data[nInf];
            nInf++;
            if (nInf + 1 == kfirst) {
              while (nPk < qEnd) {
                c_k++;
                iPk_data[c_k] = sortIdx_data[nPk - 1];
                nPk++;
              }
            }
          } else {
            iPk_data[c_k] = i1;
            nPk++;
            if (nPk == qEnd) {
              while (nInf + 1 < kfirst) {
                c_k++;
                iPk_data[c_k] = sortIdx_data[nInf];
                nInf++;
              }
            }
          }
          c_k++;
        }
        for (b_k = 0; b_k < kEnd; b_k++) {
          sortIdx_data[(j + b_k) - 1] = iPk_data[b_k];
        }
        j = qEnd;
      }
      b_i = i2;
    }
    kfirst = idelete->size[0];
    idelete->size[0] = c->size[0];
    emxEnsureCapacity_boolean_T(idelete, kfirst);
    idelete_data = idelete->data;
    for (k = 0; k < b_loop_ub; k++) {
      idelete_data[k] = false;
    }
    for (k = 0; k < b_loop_ub; k++) {
      if (!idelete_data[k]) {
        for (b_k = 0; b_k < b_loop_ub; b_k++) {
          unsigned int b_idelete_tmp;
          unsigned int idelete_tmp;
          kfirst = (b_loop_ub - b_k) - 1;
          idelete_tmp = x_data[c_data[sortIdx_data[kfirst] - 1] - 1];
          b_idelete_tmp = x_data[c_data[sortIdx_data[k] - 1] - 1];
          idelete_data[kfirst] =
              (idelete_data[kfirst] ||
               ((idelete_tmp >= (double)b_idelete_tmp - varargin_2) &&
                (idelete_tmp <= (double)b_idelete_tmp + varargin_2)));
        }
        idelete_data[k] = false;
      }
    }
    nInf = 0;
    for (k = 0; k < b_loop_ub; k++) {
      if (!idelete_data[k]) {
        nInf++;
      }
    }
    kfirst = r->size[0];
    r->size[0] = nInf;
    emxEnsureCapacity_int32_T(r, kfirst);
    iPk_data = r->data;
    kfirst = 0;
    for (k = 0; k < b_loop_ub; k++) {
      if (!idelete_data[k]) {
        iPk_data[kfirst] = k;
        kfirst++;
      }
    }
    nInf = r->size[0];
    kfirst = idx->size[0];
    idx->size[0] = r->size[0];
    emxEnsureCapacity_int32_T(idx, kfirst);
    idx_data = idx->data;
    for (k = 0; k < nInf; k++) {
      idx_data[k] = sortIdx_data[iPk_data[k]];
    }
    if (r->size[0] != 0) {
      introsort(idx, r->size[0]);
      idx_data = idx->data;
    }
  }
  emxFree_int32_T(&iwork);
  emxFree_int32_T(&y);
  emxFree_int32_T(&r);
  emxFree_boolean_T(&idelete);
  emxFree_uint32_T(&x);
  if (idx->size[0] > Yin->size[0]) {
    kfirst = sortIdx->size[0];
    sortIdx->size[0] = loop_ub;
    emxEnsureCapacity_int32_T(sortIdx, kfirst);
    kfirst = idx->size[0];
    idx->size[0] = loop_ub;
    emxEnsureCapacity_int32_T(idx, kfirst);
    idx_data = idx->data;
  } else {
    kfirst = sortIdx->size[0];
    sortIdx->size[0] = idx->size[0];
    emxEnsureCapacity_int32_T(sortIdx, kfirst);
  }
  nInf = sortIdx->size[0];
  kfirst = sortIdx->size[0];
  sortIdx->size[0] = nInf;
  emxEnsureCapacity_int32_T(sortIdx, kfirst);
  sortIdx_data = sortIdx->data;
  for (k = 0; k < nInf; k++) {
    sortIdx_data[k] = c_data[idx_data[k] - 1];
  }
  emxFree_int32_T(&c);
  emxFree_int32_T(&idx);
  kfirst = Ypk->size[0];
  Ypk->size[0] = nInf;
  emxEnsureCapacity_real_T(Ypk, kfirst);
  Ypk_data = Ypk->data;
  kfirst = Xpk->size[0];
  Xpk->size[0] = nInf;
  emxEnsureCapacity_real_T(Xpk, kfirst);
  Xpk_data = Xpk->data;
  for (k = 0; k < nInf; k++) {
    Ypk_data[k] = Yin_data[sortIdx_data[k] - 1];
    Xpk_data[k] = (unsigned int)sortIdx_data[k];
  }
  emxFree_int32_T(&sortIdx);
}

/* End of code generation (findpeaks.c) */
