//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// eml_setop.cpp
//
// Code generation for function 'eml_setop'
//

// Include files
#include "eml_setop.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
namespace coder {
void do_vectors(const array<int, 1U> &a, const array<int, 1U> &b,
                array<int, 1U> &c, array<int, 1U> &ia, array<int, 1U> &ib)
{
  int b_ialast;
  int iafirst;
  int ialast;
  int ibfirst;
  int iblast;
  int na;
  int nb;
  int nc;
  int ncmax;
  int nia;
  int nib;
  na = a.size(0);
  nb = b.size(0);
  ncmax = a.size(0) + b.size(0);
  c.set_size(ncmax);
  ia.set_size(a.size(0));
  ib.set_size(b.size(0));
  nc = -1;
  nia = -1;
  nib = 0;
  iafirst = 1;
  ialast = 0;
  ibfirst = 0;
  iblast = 0;
  while ((ialast + 1 <= na) && (iblast + 1 <= nb)) {
    int ak;
    int b_iblast;
    int bk;
    b_ialast = ialast + 1;
    ak = a[ialast];
    while ((b_ialast < a.size(0)) && (a[b_ialast] == ak)) {
      b_ialast++;
    }
    ialast = b_ialast - 1;
    b_iblast = iblast + 1;
    bk = b[iblast];
    while ((b_iblast < b.size(0)) && (b[b_iblast] == bk)) {
      b_iblast++;
    }
    iblast = b_iblast - 1;
    if (ak == bk) {
      nc++;
      c[nc] = ak;
      nia++;
      ia[nia] = iafirst;
      ialast = b_ialast;
      iafirst = b_ialast + 1;
      iblast = b_iblast;
      ibfirst = b_iblast;
    } else if (ak < bk) {
      nc++;
      nia++;
      c[nc] = ak;
      ia[nia] = iafirst;
      ialast = b_ialast;
      iafirst = b_ialast + 1;
    } else {
      nc++;
      nib++;
      c[nc] = bk;
      ib[nib - 1] = ibfirst + 1;
      iblast = b_iblast;
      ibfirst = b_iblast;
    }
  }
  while (ialast + 1 <= na) {
    b_ialast = ialast + 1;
    while ((b_ialast < a.size(0)) && (a[b_ialast] == a[ialast])) {
      b_ialast++;
    }
    nc++;
    nia++;
    c[nc] = a[ialast];
    ia[nia] = iafirst;
    ialast = b_ialast;
    iafirst = b_ialast + 1;
  }
  while (iblast + 1 <= nb) {
    b_ialast = iblast + 1;
    while ((b_ialast < b.size(0)) && (b[b_ialast] == b[iblast])) {
      b_ialast++;
    }
    nc++;
    nib++;
    c[nc] = b[iblast];
    ib[nib - 1] = ibfirst + 1;
    iblast = b_ialast;
    ibfirst = b_ialast;
  }
  if (a.size(0) > 0) {
    if (nia + 1 < 1) {
      b_ialast = 0;
    } else {
      b_ialast = nia + 1;
    }
    ia.set_size(b_ialast);
  }
  if (b.size(0) > 0) {
    if (nib < 1) {
      nib = 0;
    }
    ib.set_size(nib);
  }
  if (ncmax > 0) {
    if (nc + 1 < 1) {
      b_ialast = 0;
    } else {
      b_ialast = nc + 1;
    }
    c.set_size(b_ialast);
  }
}

} // namespace coder

// End of code generation (eml_setop.cpp)
