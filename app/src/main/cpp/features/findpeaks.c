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
#include "introsort.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
int findpeaks(const double Yin_data[], int Yin_size, double varargin_2,
              double Ypk_data[], double Xpk_data[], int *Xpk_size)
{
  double ykfirst;
  int c_data[300];
  int iwork_data[300];
  int fPk_data[150];
  int iFinite_data[150];
  int iInfinite_data[150];
  int iPk_data[150];
  int Ypk_size;
  int b_nPk;
  int jj;
  int k;
  int kfirst;
  int nInf;
  int nPk;
  int qEnd;
  short y_data[300];
  unsigned char x_data[150];
  char dir;
  boolean_T idelete_data[300];
  boolean_T isinfykfirst;
  kfirst = Yin_size - 1;
  for (k = 0; k <= kfirst; k++) {
    x_data[k] = (unsigned char)((unsigned int)k + 1U);
  }
  nPk = 0;
  nInf = 0;
  dir = 'n';
  kfirst = 0;
  ykfirst = rtInf;
  isinfykfirst = true;
  for (k = 1; k <= Yin_size; k++) {
    double yk;
    boolean_T isinfyk;
    yk = Yin_data[k - 1];
    if (rtIsNaN(yk)) {
      yk = rtInf;
      isinfyk = true;
    } else if (rtIsInf(yk) && (yk > 0.0)) {
      isinfyk = true;
      nInf++;
      iInfinite_data[nInf - 1] = k;
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
          iFinite_data[nPk - 1] = kfirst;
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
  b_nPk = 0;
  for (k = 0; k < nPk; k++) {
    kfirst = iFinite_data[k];
    ykfirst = Yin_data[kfirst - 1];
    if ((ykfirst > rtMinusInf) &&
        (ykfirst - fmax(Yin_data[kfirst - 2], Yin_data[kfirst]) >= 0.0)) {
      b_nPk++;
      iPk_data[b_nPk - 1] = kfirst;
    }
  }
  if (nInf < 1) {
    nInf = 0;
  }
  Ypk_size = do_vectors(iPk_data, b_nPk, iInfinite_data, nInf, c_data, fPk_data,
                        &Ypk_size, iFinite_data, &kfirst);
  if ((Ypk_size == 0) || (varargin_2 == 0.0)) {
    if (Ypk_size > 0) {
      y_data[0] = 1;
      kfirst = 1;
      for (k = 2; k <= Ypk_size; k++) {
        kfirst++;
        y_data[k - 1] = (short)kfirst;
      }
    }
    for (k = 0; k < Ypk_size; k++) {
      iwork_data[k] = y_data[k];
    }
  } else {
    int sortIdx_data[300];
    int i;
    kfirst = Ypk_size - 1;
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
    if (((unsigned int)Ypk_size & 1U) != 0U) {
      sortIdx_data[Ypk_size - 1] = Ypk_size;
    }
    i = 2;
    while (i < Ypk_size) {
      int i2;
      int j;
      i2 = i << 1;
      j = 1;
      for (kfirst = i + 1; kfirst < Ypk_size + 1; kfirst = qEnd + i) {
        int kEnd;
        nPk = j - 1;
        b_nPk = kfirst;
        qEnd = j + i2;
        if (qEnd > Ypk_size + 1) {
          qEnd = Ypk_size + 1;
        }
        nInf = 0;
        kEnd = qEnd - j;
        while (nInf < kEnd) {
          int b_i;
          ykfirst = Yin_data[c_data[sortIdx_data[nPk] - 1] - 1];
          b_i = sortIdx_data[b_nPk - 1];
          if ((ykfirst >= Yin_data[c_data[b_i - 1] - 1]) || rtIsNaN(ykfirst)) {
            iwork_data[nInf] = sortIdx_data[nPk];
            nPk++;
            if (nPk + 1 == kfirst) {
              while (b_nPk < qEnd) {
                nInf++;
                iwork_data[nInf] = sortIdx_data[b_nPk - 1];
                b_nPk++;
              }
            }
          } else {
            iwork_data[nInf] = b_i;
            b_nPk++;
            if (b_nPk == qEnd) {
              while (nPk + 1 < kfirst) {
                nInf++;
                iwork_data[nInf] = sortIdx_data[nPk];
                nPk++;
              }
            }
          }
          nInf++;
        }
        for (k = 0; k < kEnd; k++) {
          sortIdx_data[(j + k) - 1] = iwork_data[k];
        }
        j = qEnd;
      }
      i = i2;
    }
    memset(&idelete_data[0], 0, (unsigned int)Ypk_size * sizeof(boolean_T));
    for (k = 0; k < Ypk_size; k++) {
      if (!idelete_data[k]) {
        kfirst = x_data[c_data[sortIdx_data[k] - 1] - 1];
        for (jj = 0; jj < Ypk_size; jj++) {
          nPk = (Ypk_size - jj) - 1;
          b_nPk = x_data[c_data[sortIdx_data[nPk] - 1] - 1];
          idelete_data[nPk] =
              (idelete_data[nPk] || ((b_nPk >= (double)kfirst - varargin_2) &&
                                     (b_nPk <= (double)kfirst + varargin_2)));
        }
        idelete_data[k] = false;
      }
    }
    nPk = 0;
    kfirst = 0;
    for (k = 0; k < Ypk_size; k++) {
      if (!idelete_data[k]) {
        nPk++;
        y_data[kfirst] = (short)k;
        kfirst++;
      }
    }
    Ypk_size = nPk;
    for (k = 0; k < nPk; k++) {
      iwork_data[k] = sortIdx_data[y_data[k]];
    }
    if (nPk != 0) {
      introsort(iwork_data, nPk);
    }
  }
  if (Ypk_size > Yin_size) {
    Ypk_size = Yin_size;
  }
  *Xpk_size = Ypk_size;
  for (k = 0; k < Ypk_size; k++) {
    kfirst = c_data[iwork_data[k] - 1];
    Ypk_data[k] = Yin_data[kfirst - 1];
    Xpk_data[k] = (unsigned char)((unsigned char)(kfirst - 1) + 1);
  }
  return Ypk_size;
}

/* End of code generation (findpeaks.c) */
