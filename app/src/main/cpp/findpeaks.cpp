//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// findpeaks.cpp
//
// Code generation for function 'findpeaks'
//

// Include files
#include "findpeaks.h"
#include "eml_setop.h"
#include "introsort.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>
#include <cstring>

// Function Definitions
namespace coder {
void findpeaks(const array<double, 1U> &Yin, array<double, 1U> &Ypk,
               array<double, 1U> &Xpk)
{
  array<int, 1U> c;
  array<int, 1U> iInfinite;
  array<int, 1U> iPk;
  array<int, 1U> idx;
  array<int, 1U> iwork;
  array<int, 1U> r;
  array<int, 1U> sortIdx;
  array<unsigned int, 1U> x;
  array<bool, 1U> idelete;
  double ykfirst;
  int kfirst;
  int nInf;
  int nPk;
  int q;
  int qEnd;
  char dir;
  bool isinfykfirst;
  q = Yin.size(0);
  x.set_size(Yin.size(0));
  for (int k{0}; k < q; k++) {
    x[k] = static_cast<unsigned int>(k) + 1U;
  }
  idx.set_size(Yin.size(0));
  iInfinite.set_size(Yin.size(0));
  nPk = 0;
  nInf = 0;
  dir = 'n';
  kfirst = 0;
  ykfirst = rtInf;
  isinfykfirst = true;
  for (int k{1}; k <= q; k++) {
    double yk;
    bool isinfyk;
    yk = Yin[k - 1];
    if (std::isnan(yk)) {
      yk = rtInf;
      isinfyk = true;
    } else if (std::isinf(yk) && (yk > 0.0)) {
      isinfyk = true;
      nInf++;
      iInfinite[nInf - 1] = k;
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
          idx[nPk - 1] = kfirst;
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
  idx.set_size(nPk);
  if (nInf < 1) {
    nInf = 0;
  }
  iInfinite.set_size(nInf);
  iPk.set_size(nPk);
  nInf = 0;
  for (int k{0}; k < nPk; k++) {
    ykfirst = Yin[idx[k] - 1];
    if ((ykfirst > rtMinusInf) &&
        (ykfirst - std::fmax(Yin[idx[k] - 2], Yin[idx[k]]) >= 0.0)) {
      nInf++;
      iPk[nInf - 1] = idx[k];
    }
  }
  iPk.set_size(nInf);
  do_vectors(iPk, iInfinite, c, sortIdx, idx);
  if (c.size(0) == 0) {
    idx.set_size(0);
  } else {
    int i;
    int loop_ub;
    int n;
    loop_ub = c.size(0);
    sortIdx.set_size(c.size(0));
    n = c.size(0) + 1;
    iwork.set_size(c.size(0));
    nInf = c.size(0) - 1;
    for (int k{1}; k <= nInf; k += 2) {
      ykfirst = Yin[c[k - 1] - 1];
      if ((ykfirst >= Yin[c[k] - 1]) || std::isnan(ykfirst)) {
        sortIdx[k - 1] = k;
        sortIdx[k] = k + 1;
      } else {
        sortIdx[k - 1] = k + 1;
        sortIdx[k] = k;
      }
    }
    if ((static_cast<unsigned int>(c.size(0)) & 1U) != 0U) {
      sortIdx[c.size(0) - 1] = c.size(0);
    }
    i = 2;
    while (i < n - 1) {
      int i2;
      int j;
      i2 = i << 1;
      j = 1;
      for (nInf = i + 1; nInf < n; nInf = qEnd + i) {
        int kEnd;
        kfirst = j - 1;
        q = nInf;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }
        nPk = 0;
        kEnd = qEnd - j;
        while (nPk < kEnd) {
          int b_i;
          ykfirst = Yin[c[sortIdx[kfirst] - 1] - 1];
          b_i = sortIdx[q - 1];
          if ((ykfirst >= Yin[c[b_i - 1] - 1]) || std::isnan(ykfirst)) {
            iwork[nPk] = sortIdx[kfirst];
            kfirst++;
            if (kfirst + 1 == nInf) {
              while (q < qEnd) {
                nPk++;
                iwork[nPk] = sortIdx[q - 1];
                q++;
              }
            }
          } else {
            iwork[nPk] = b_i;
            q++;
            if (q == qEnd) {
              while (kfirst + 1 < nInf) {
                nPk++;
                iwork[nPk] = sortIdx[kfirst];
                kfirst++;
              }
            }
          }
          nPk++;
        }
        for (int b_k{0}; b_k < kEnd; b_k++) {
          sortIdx[(j + b_k) - 1] = iwork[b_k];
        }
        j = qEnd;
      }
      i = i2;
    }
    idelete.set_size(c.size(0));
    std::memset(&idelete[0], 0,
                static_cast<unsigned int>(loop_ub) * sizeof(bool));
    for (int k{0}; k < loop_ub; k++) {
      if (!idelete[k]) {
        unsigned int idelete_tmp;
        idelete_tmp = x[c[sortIdx[k] - 1] - 1];
        for (int b_k{0}; b_k < loop_ub; b_k++) {
          unsigned int b_idelete_tmp;
          nInf = (loop_ub - b_k) - 1;
          b_idelete_tmp = x[c[sortIdx[nInf] - 1] - 1];
          idelete[nInf] =
              (idelete[nInf] || ((static_cast<int>(b_idelete_tmp) >=
                                  static_cast<int>(idelete_tmp) - 25) &&
                                 (b_idelete_tmp <= idelete_tmp + 25U)));
        }
        idelete[k] = false;
      }
    }
    nInf = 0;
    for (int k{0}; k < loop_ub; k++) {
      if (!idelete[k]) {
        nInf++;
      }
    }
    r.set_size(nInf);
    nInf = 0;
    for (int k{0}; k < loop_ub; k++) {
      if (!idelete[k]) {
        r[nInf] = k;
        nInf++;
      }
    }
    nInf = r.size(0);
    idx.set_size(r.size(0));
    for (int k{0}; k < nInf; k++) {
      idx[k] = sortIdx[r[k]];
    }
    if (r.size(0) != 0) {
      internal::introsort(idx, r.size(0));
    }
  }
  if (idx.size(0) > Yin.size(0)) {
    sortIdx.set_size(Yin.size(0));
    idx.set_size(Yin.size(0));
  } else {
    sortIdx.set_size(idx.size(0));
  }
  nInf = sortIdx.size(0);
  sortIdx.set_size(sortIdx.size(0));
  for (int k{0}; k < nInf; k++) {
    sortIdx[k] = c[idx[k] - 1];
  }
  Ypk.set_size(nInf);
  Xpk.set_size(nInf);
  for (int k{0}; k < nInf; k++) {
    Ypk[k] = Yin[sortIdx[k] - 1];
    Xpk[k] = static_cast<unsigned int>(sortIdx[k]);
  }
}

} // namespace coder

// End of code generation (findpeaks.cpp)
