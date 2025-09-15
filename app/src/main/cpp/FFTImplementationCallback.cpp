//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// FFTImplementationCallback.cpp
//
// Code generation for function 'FFTImplementationCallback'
//

// Include files
#include "FFTImplementationCallback.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>
#include <cstring>

// Function Definitions
namespace coder {
namespace internal {
namespace fft {
void FFTImplementationCallback::doHalfLengthBluestein(
    const array<double, 1U> &x, array<creal_T, 1U> &y, int nrowsx, int nRows,
    int nfft, const array<creal_T, 1U> &wwc, const array<double, 2U> &costab,
    const array<double, 2U> &sintab, const array<double, 2U> &costabinv,
    const array<double, 2U> &sintabinv)
{
  array<creal_T, 1U> fv;
  array<creal_T, 1U> fy;
  array<creal_T, 1U> reconVar1;
  array<creal_T, 1U> reconVar2;
  array<creal_T, 1U> ytmp;
  array<double, 2U> b_costab;
  array<double, 2U> b_sintab;
  array<double, 2U> costab1q;
  array<double, 2U> hcostabinv;
  array<double, 2U> hsintab;
  array<double, 2U> hsintabinv;
  array<int, 2U> wrapIndex;
  double e;
  double im;
  double temp_im;
  double temp_re;
  double twid_im;
  double twid_re;
  int b_k;
  int b_nfft;
  int hnRows;
  int iDelta;
  int iDelta2;
  int iheight;
  int ju;
  int n;
  int n2;
  int nRowsD2;
  int nd2;
  bool tst;
  hnRows = static_cast<int>(static_cast<unsigned int>(nRows) >> 1);
  ytmp.set_size(hnRows);
  if (hnRows > nrowsx) {
    ytmp.set_size(hnRows);
    if (hnRows - 1 >= 0) {
      std::memset(&ytmp[0], 0,
                  static_cast<unsigned int>(hnRows) * sizeof(creal_T));
    }
  }
  if ((static_cast<unsigned int>(x.size(0)) & 1U) == 0U) {
    tst = true;
    ju = x.size(0);
  } else if (x.size(0) >= nRows) {
    tst = true;
    ju = nRows;
  } else {
    tst = false;
    ju = x.size(0) - 1;
  }
  nd2 = nRows << 1;
  e = 6.2831853071795862 / static_cast<double>(nd2);
  n2 = static_cast<int>((static_cast<unsigned int>(nd2) >> 1) >> 1);
  costab1q.set_size(1, n2 + 1);
  costab1q[0] = 1.0;
  nd2 = static_cast<int>(static_cast<unsigned int>(n2) >> 1) - 1;
  for (int k{0}; k <= nd2; k++) {
    costab1q[k + 1] = std::cos(e * (static_cast<double>(k) + 1.0));
  }
  n = nd2 + 2;
  for (int k{n}; k < n2; k++) {
    costab1q[k] = std::sin(e * static_cast<double>(n2 - k));
  }
  costab1q[n2] = 0.0;
  n = costab1q.size(1) - 1;
  n2 = (costab1q.size(1) - 1) << 1;
  b_costab.set_size(1, n2 + 1);
  b_sintab.set_size(1, n2 + 1);
  b_costab[0] = 1.0;
  b_sintab[0] = 0.0;
  for (int k{0}; k < n; k++) {
    b_costab[k + 1] = costab1q[k + 1];
    b_sintab[k + 1] = -costab1q[(n - k) - 1];
  }
  nd2 = costab1q.size(1);
  for (int k{nd2}; k <= n2; k++) {
    b_costab[k] = -costab1q[n2 - k];
    b_sintab[k] = -costab1q[k - n];
  }
  nd2 = static_cast<int>(static_cast<unsigned int>(costab.size(1)) >> 1);
  costab1q.set_size(1, nd2);
  hsintab.set_size(1, nd2);
  hcostabinv.set_size(1, nd2);
  hsintabinv.set_size(1, nd2);
  for (int k{0}; k < nd2; k++) {
    n = ((k + 1) << 1) - 2;
    costab1q[k] = costab[n];
    hsintab[k] = sintab[n];
    hcostabinv[k] = costabinv[n];
    hsintabinv[k] = sintabinv[n];
  }
  reconVar1.set_size(hnRows);
  reconVar2.set_size(hnRows);
  wrapIndex.set_size(1, hnRows);
  for (int k{0}; k < hnRows; k++) {
    nd2 = k << 1;
    e = b_sintab[nd2];
    im = b_costab[nd2];
    reconVar1[k].re = e + 1.0;
    reconVar1[k].im = -im;
    reconVar2[k].re = 1.0 - e;
    reconVar2[k].im = im;
    if (k != 0) {
      wrapIndex[k] = (hnRows - k) + 1;
    } else {
      wrapIndex[0] = 1;
    }
  }
  if (ju > nRows) {
    ju = nRows;
  }
  n = static_cast<int>(static_cast<unsigned int>(ju) >> 1) - 1;
  for (int k{0}; k <= n; k++) {
    nd2 = (hnRows + k) - 1;
    temp_re = wwc[nd2].re;
    temp_im = wwc[nd2].im;
    nd2 = k << 1;
    twid_re = x[nd2];
    twid_im = x[nd2 + 1];
    ytmp[k].re = temp_re * twid_re + temp_im * twid_im;
    ytmp[k].im = temp_re * twid_im - temp_im * twid_re;
  }
  if (!tst) {
    nd2 = hnRows + n;
    temp_re = wwc[nd2].re;
    temp_im = wwc[nd2].im;
    if (n < 0) {
      nd2 = 0;
    } else {
      nd2 = (n + 1) << 1;
    }
    twid_re = x[nd2];
    ytmp[n + 1].re = temp_re * twid_re + temp_im * 0.0;
    ytmp[n + 1].im = temp_re * 0.0 - temp_im * twid_re;
    if (n + 3 <= hnRows) {
      nd2 = n + 3;
      if (nd2 <= hnRows) {
        std::memset(&ytmp[nd2 + -1], 0,
                    static_cast<unsigned int>((hnRows - nd2) + 1) *
                        sizeof(creal_T));
      }
    }
  } else if (n + 2 <= hnRows) {
    nd2 = n + 2;
    if (nd2 <= hnRows) {
      std::memset(&ytmp[nd2 + -1], 0,
                  static_cast<unsigned int>((hnRows - nd2) + 1) *
                      sizeof(creal_T));
    }
  }
  b_nfft = static_cast<int>(static_cast<unsigned int>(nfft) >> 1);
  fy.set_size(b_nfft);
  if (b_nfft > ytmp.size(0)) {
    fy.set_size(b_nfft);
    std::memset(&fy[0], 0, static_cast<unsigned int>(b_nfft) * sizeof(creal_T));
  }
  n = ytmp.size(0);
  if (n > b_nfft) {
    n = b_nfft;
  }
  iDelta = b_nfft - 2;
  nRowsD2 = static_cast<int>(static_cast<unsigned int>(b_nfft) >> 1);
  b_k = static_cast<int>(static_cast<unsigned int>(nRowsD2) >> 1);
  n2 = 0;
  ju = 0;
  for (int k{0}; k <= n - 2; k++) {
    fy[n2] = ytmp[k];
    nd2 = b_nfft;
    tst = true;
    while (tst) {
      nd2 >>= 1;
      ju = static_cast<int>(static_cast<unsigned int>(ju) ^
                            static_cast<unsigned int>(nd2));
      tst = ((static_cast<unsigned int>(ju) & static_cast<unsigned int>(nd2)) ==
             0U);
    }
    n2 = ju;
  }
  if (n - 2 < 0) {
    nd2 = 0;
  } else {
    nd2 = n - 1;
  }
  fy[n2] = ytmp[nd2];
  if (b_nfft > 1) {
    for (int k{0}; k <= iDelta; k += 2) {
      temp_re = fy[k + 1].re;
      temp_im = fy[k + 1].im;
      e = fy[k].re;
      im = fy[k].im;
      fy[k + 1].re = e - temp_re;
      fy[k + 1].im = fy[k].im - fy[k + 1].im;
      e += temp_re;
      im += temp_im;
      fy[k].re = e;
      fy[k].im = im;
    }
  }
  iDelta = 2;
  iDelta2 = 4;
  iheight = ((b_k - 1) << 2) + 1;
  while (b_k > 0) {
    int i;
    for (i = 0; i < iheight; i += iDelta2) {
      nd2 = i + iDelta;
      temp_re = fy[nd2].re;
      temp_im = fy[nd2].im;
      fy[nd2].re = fy[i].re - temp_re;
      fy[nd2].im = fy[i].im - temp_im;
      fy[i].re = fy[i].re + temp_re;
      fy[i].im = fy[i].im + temp_im;
    }
    nd2 = 1;
    for (n = b_k; n < nRowsD2; n += b_k) {
      twid_re = costab1q[n];
      twid_im = hsintab[n];
      i = nd2;
      n2 = nd2 + iheight;
      while (i < n2) {
        ju = i + iDelta;
        e = fy[ju].im;
        im = fy[ju].re;
        temp_re = twid_re * im - twid_im * e;
        temp_im = twid_re * e + twid_im * im;
        fy[ju].re = fy[i].re - temp_re;
        fy[ju].im = fy[i].im - temp_im;
        fy[i].re = fy[i].re + temp_re;
        fy[i].im = fy[i].im + temp_im;
        i += iDelta2;
      }
      nd2++;
    }
    b_k = static_cast<int>(static_cast<unsigned int>(b_k) >> 1);
    iDelta = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iDelta;
  }
  FFTImplementationCallback::r2br_r2dit_trig_impl(wwc, b_nfft, costab1q,
                                                  hsintab, fv);
  for (int k{0}; k < b_nfft; k++) {
    e = fy[k].re;
    im = fv[k].im;
    temp_im = fy[k].im;
    twid_re = fv[k].re;
    fy[k].re = e * twid_re - temp_im * im;
    fy[k].im = e * im + temp_im * twid_re;
  }
  FFTImplementationCallback::r2br_r2dit_trig_impl(fy, b_nfft, hcostabinv,
                                                  hsintabinv, fv);
  if (fv.size(0) > 1) {
    e = 1.0 / static_cast<double>(fv.size(0));
    nd2 = fv.size(0);
    for (int k{0}; k < nd2; k++) {
      fv[k].re = e * fv[k].re;
      fv[k].im = e * fv[k].im;
    }
  }
  nd2 = wwc.size(0);
  for (int k{hnRows}; k <= nd2; k++) {
    e = wwc[k - 1].re;
    im = fv[k - 1].im;
    temp_im = wwc[k - 1].im;
    twid_re = fv[k - 1].re;
    n = k - hnRows;
    ytmp[n].re = e * twid_re + temp_im * im;
    ytmp[n].im = e * im - temp_im * twid_re;
  }
  for (int k{0}; k < hnRows; k++) {
    double b_ytmp_re_tmp;
    double ytmp_re_tmp;
    nd2 = wrapIndex[k];
    e = ytmp[k].re;
    im = reconVar1[k].im;
    temp_im = ytmp[k].im;
    twid_re = reconVar1[k].re;
    twid_im = ytmp[nd2 - 1].re;
    temp_re = -ytmp[nd2 - 1].im;
    ytmp_re_tmp = reconVar2[k].im;
    b_ytmp_re_tmp = reconVar2[k].re;
    y[k].re = 0.5 * ((e * twid_re - temp_im * im) +
                     (twid_im * b_ytmp_re_tmp - temp_re * ytmp_re_tmp));
    y[k].im = 0.5 * ((e * im + temp_im * twid_re) +
                     (twid_im * ytmp_re_tmp + temp_re * b_ytmp_re_tmp));
    nd2 = hnRows + k;
    y[nd2].re = 0.5 * ((e * b_ytmp_re_tmp - temp_im * ytmp_re_tmp) +
                       (twid_im * twid_re - temp_re * im));
    y[nd2].im = 0.5 * ((e * ytmp_re_tmp + temp_im * b_ytmp_re_tmp) +
                       (twid_im * im + temp_re * twid_re));
  }
}

void FFTImplementationCallback::r2br_r2dit_trig_impl(
    const array<creal_T, 1U> &x, int unsigned_nRows,
    const array<double, 2U> &costab, const array<double, 2U> &sintab,
    array<creal_T, 1U> &y)
{
  double im;
  double re;
  double temp_im;
  double temp_re;
  int iDelta;
  int iDelta2;
  int iheight;
  int istart;
  int iy;
  int j;
  int ju;
  int k;
  int nRowsD2;
  y.set_size(unsigned_nRows);
  if (unsigned_nRows > x.size(0)) {
    y.set_size(unsigned_nRows);
    for (int i{0}; i < unsigned_nRows; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  }
  j = x.size(0);
  if (j > unsigned_nRows) {
    j = unsigned_nRows;
  }
  iDelta = unsigned_nRows - 2;
  nRowsD2 = static_cast<int>(static_cast<unsigned int>(unsigned_nRows) >> 1);
  k = static_cast<int>(static_cast<unsigned int>(nRowsD2) >> 1);
  iy = 0;
  ju = 0;
  for (int i{0}; i <= j - 2; i++) {
    bool tst;
    y[iy] = x[i];
    istart = unsigned_nRows;
    tst = true;
    while (tst) {
      istart >>= 1;
      ju ^= istart;
      tst = ((ju & istart) == 0);
    }
    iy = ju;
  }
  if (j - 2 < 0) {
    istart = 0;
  } else {
    istart = j - 1;
  }
  y[iy] = x[istart];
  if (unsigned_nRows > 1) {
    for (int i{0}; i <= iDelta; i += 2) {
      temp_re = y[i + 1].re;
      temp_im = y[i + 1].im;
      re = y[i].re;
      im = y[i].im;
      y[i + 1].re = re - temp_re;
      y[i + 1].im = y[i].im - y[i + 1].im;
      re += temp_re;
      im += temp_im;
      y[i].re = re;
      y[i].im = im;
    }
  }
  iDelta = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int b_i;
    for (b_i = 0; b_i < iheight; b_i += iDelta2) {
      istart = b_i + iDelta;
      temp_re = y[istart].re;
      temp_im = y[istart].im;
      y[istart].re = y[b_i].re - temp_re;
      y[istart].im = y[b_i].im - temp_im;
      y[b_i].re = y[b_i].re + temp_re;
      y[b_i].im = y[b_i].im + temp_im;
    }
    istart = 1;
    for (j = k; j < nRowsD2; j += k) {
      double twid_im;
      double twid_re;
      twid_re = costab[j];
      twid_im = sintab[j];
      b_i = istart;
      iy = istart + iheight;
      while (b_i < iy) {
        ju = b_i + iDelta;
        re = y[ju].im;
        im = y[ju].re;
        temp_re = twid_re * im - twid_im * re;
        temp_im = twid_re * re + twid_im * im;
        y[ju].re = y[b_i].re - temp_re;
        y[ju].im = y[b_i].im - temp_im;
        y[b_i].re = y[b_i].re + temp_re;
        y[b_i].im = y[b_i].im + temp_im;
        b_i += iDelta2;
      }
      istart++;
    }
    k = static_cast<int>(static_cast<unsigned int>(k) >> 1);
    iDelta = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iDelta;
  }
}

void FFTImplementationCallback::doHalfLengthRadix2(
    const array<double, 1U> &x, array<creal_T, 1U> &y, int unsigned_nRows,
    const array<double, 2U> &costab, const array<double, 2U> &sintab)
{
  array<creal_T, 1U> reconVar1;
  array<creal_T, 1U> reconVar2;
  array<double, 2U> hcostab;
  array<double, 2U> hsintab;
  array<int, 2U> wrapIndex;
  array<int, 1U> bitrevIndex;
  double b_im;
  double b_re;
  double im;
  double re;
  double temp2_im;
  double temp2_re;
  double temp_im;
  double temp_re;
  int hszCostab;
  int iDelta;
  int iDelta2;
  int iheight;
  int ihi;
  int iy;
  int k;
  int nRows;
  int nRowsD2;
  int nRowsM1;
  bool tst;
  nRows = static_cast<int>(static_cast<unsigned int>(unsigned_nRows) >> 1);
  ihi = y.size(0);
  if (ihi > nRows) {
    ihi = nRows;
  }
  nRowsM1 = ihi - 2;
  iDelta = nRows - 2;
  nRowsD2 = static_cast<int>(static_cast<unsigned int>(nRows) >> 1);
  k = static_cast<int>(static_cast<unsigned int>(nRowsD2) >> 1);
  hszCostab = static_cast<int>(static_cast<unsigned int>(costab.size(1)) >> 1);
  hcostab.set_size(1, hszCostab);
  hsintab.set_size(1, hszCostab);
  for (int i{0}; i < hszCostab; i++) {
    iy = ((i + 1) << 1) - 2;
    hcostab[i] = costab[iy];
    hsintab[i] = sintab[iy];
  }
  reconVar1.set_size(nRows);
  reconVar2.set_size(nRows);
  wrapIndex.set_size(1, nRows);
  hszCostab = 0;
  iy = 1;
  bitrevIndex.set_size(nRows);
  for (int i{0}; i < nRows; i++) {
    re = sintab[i];
    im = costab[i];
    reconVar1[i].re = re + 1.0;
    reconVar1[i].im = -im;
    reconVar2[i].re = 1.0 - re;
    reconVar2[i].im = im;
    if (i != 0) {
      wrapIndex[i] = (nRows - i) + 1;
    } else {
      wrapIndex[0] = 1;
    }
    bitrevIndex[i] = 0;
  }
  for (int i{0}; i <= nRowsM1; i++) {
    bitrevIndex[i] = iy;
    iy = nRows;
    tst = true;
    while (tst) {
      iy >>= 1;
      hszCostab = static_cast<int>(static_cast<unsigned int>(hszCostab) ^
                                   static_cast<unsigned int>(iy));
      tst = ((static_cast<unsigned int>(hszCostab) &
              static_cast<unsigned int>(iy)) == 0U);
    }
    iy = hszCostab + 1;
  }
  bitrevIndex[ihi - 1] = iy;
  if ((static_cast<unsigned int>(x.size(0)) & 1U) == 0U) {
    tst = true;
    hszCostab = x.size(0);
  } else if (x.size(0) >= unsigned_nRows) {
    tst = true;
    hszCostab = unsigned_nRows;
  } else {
    tst = false;
    hszCostab = x.size(0) - 1;
  }
  if (hszCostab > unsigned_nRows) {
    hszCostab = unsigned_nRows;
  }
  iy = static_cast<int>(static_cast<unsigned int>(hszCostab) >> 1);
  for (int i{0}; i < iy; i++) {
    hszCostab = i << 1;
    y[bitrevIndex[i] - 1].re = x[hszCostab];
    y[bitrevIndex[i] - 1].im = x[hszCostab + 1];
  }
  if (!tst) {
    if (iy - 1 < 0) {
      hszCostab = 0;
    } else {
      hszCostab = iy << 1;
    }
    y[bitrevIndex[iy] - 1].re = x[hszCostab];
    y[bitrevIndex[iy] - 1].im = 0.0;
  }
  if (nRows > 1) {
    for (int i{0}; i <= iDelta; i += 2) {
      re = y[i + 1].re;
      im = y[i + 1].im;
      temp_re = re;
      temp_im = im;
      b_re = y[i].re;
      b_im = y[i].im;
      re = b_re - re;
      im = b_im - im;
      y[i + 1].re = re;
      y[i + 1].im = im;
      b_re += temp_re;
      b_im += temp_im;
      y[i].re = b_re;
      y[i].im = b_im;
    }
  }
  iDelta = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int b_i;
    for (b_i = 0; b_i < iheight; b_i += iDelta2) {
      hszCostab = b_i + iDelta;
      temp_re = y[hszCostab].re;
      temp_im = y[hszCostab].im;
      y[hszCostab].re = y[b_i].re - temp_re;
      y[hszCostab].im = y[b_i].im - temp_im;
      y[b_i].re = y[b_i].re + temp_re;
      y[b_i].im = y[b_i].im + temp_im;
    }
    hszCostab = 1;
    for (iy = k; iy < nRowsD2; iy += k) {
      temp2_re = hcostab[iy];
      temp2_im = hsintab[iy];
      b_i = hszCostab;
      ihi = hszCostab + iheight;
      while (b_i < ihi) {
        nRowsM1 = b_i + iDelta;
        re = y[nRowsM1].im;
        im = y[nRowsM1].re;
        temp_re = temp2_re * im - temp2_im * re;
        temp_im = temp2_re * re + temp2_im * im;
        y[nRowsM1].re = y[b_i].re - temp_re;
        y[nRowsM1].im = y[b_i].im - temp_im;
        y[b_i].re = y[b_i].re + temp_re;
        y[b_i].im = y[b_i].im + temp_im;
        b_i += iDelta2;
      }
      hszCostab++;
    }
    k = static_cast<int>(static_cast<unsigned int>(k) >> 1);
    iDelta = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iDelta;
  }
  temp_re = y[0].re;
  temp_im = y[0].im;
  im = y[0].re * reconVar1[0].re;
  b_re = y[0].re * reconVar1[0].im;
  b_im = -y[0].im;
  temp2_re = temp_re * reconVar2[0].re;
  re = temp_re * reconVar2[0].im;
  y[0].re = 0.5 * ((im - y[0].im * reconVar1[0].im) +
                   (temp2_re - b_im * reconVar2[0].im));
  y[0].im = 0.5 * ((b_re + y[0].im * reconVar1[0].re) +
                   (re + b_im * reconVar2[0].re));
  y[nRows].re = 0.5 * ((temp2_re - temp_im * reconVar2[0].im) +
                       (im - b_im * reconVar1[0].im));
  y[nRows].im = 0.5 * ((re + temp_im * reconVar2[0].re) +
                       (b_re + b_im * reconVar1[0].re));
  for (int i{2}; i <= nRowsD2; i++) {
    temp_re = y[i - 1].re;
    temp_im = y[i - 1].im;
    hszCostab = wrapIndex[i - 1];
    temp2_re = y[hszCostab - 1].re;
    temp2_im = y[hszCostab - 1].im;
    re = reconVar1[i - 1].im;
    im = reconVar1[i - 1].re;
    b_re = reconVar2[i - 1].im;
    b_im = reconVar2[i - 1].re;
    y[i - 1].re = 0.5 * ((temp_re * im - temp_im * re) +
                         (temp2_re * b_im - -temp2_im * b_re));
    y[i - 1].im = 0.5 * ((temp_re * re + temp_im * im) +
                         (temp2_re * b_re + -temp2_im * b_im));
    iy = (nRows + i) - 1;
    y[iy].re = 0.5 * ((temp_re * b_im - temp_im * b_re) +
                      (temp2_re * im - -temp2_im * re));
    y[iy].im = 0.5 * ((temp_re * b_re + temp_im * b_im) +
                      (temp2_re * re + -temp2_im * im));
    re = reconVar1[hszCostab - 1].im;
    im = reconVar1[hszCostab - 1].re;
    b_re = reconVar2[hszCostab - 1].im;
    b_im = reconVar2[hszCostab - 1].re;
    y[hszCostab - 1].re = 0.5 * ((temp2_re * im - temp2_im * re) +
                                 (temp_re * b_im - -temp_im * b_re));
    y[hszCostab - 1].im = 0.5 * ((temp2_re * re + temp2_im * im) +
                                 (temp_re * b_re + -temp_im * b_im));
    hszCostab = (hszCostab + nRows) - 1;
    y[hszCostab].re = 0.5 * ((temp2_re * b_im - temp2_im * b_re) +
                             (temp_re * im - -temp_im * re));
    y[hszCostab].im = 0.5 * ((temp2_re * b_re + temp2_im * b_im) +
                             (temp_re * re + -temp_im * im));
  }
  double b_y_re_tmp;
  double y_re_tmp;
  temp_re = y[nRowsD2].re;
  temp_im = y[nRowsD2].im;
  im = reconVar1[nRowsD2].im;
  b_re = reconVar1[nRowsD2].re;
  b_im = temp_re * b_re;
  temp2_re = temp_re * im;
  temp2_im = reconVar2[nRowsD2].im;
  y_re_tmp = reconVar2[nRowsD2].re;
  b_y_re_tmp = temp_re * y_re_tmp;
  re = temp_re * temp2_im;
  y[nRowsD2].re =
      0.5 * ((b_im - temp_im * im) + (b_y_re_tmp - -temp_im * temp2_im));
  y[nRowsD2].im =
      0.5 * ((temp2_re + temp_im * b_re) + (re + -temp_im * y_re_tmp));
  hszCostab = nRows + nRowsD2;
  y[hszCostab].re =
      0.5 * ((b_y_re_tmp - temp_im * temp2_im) + (b_im - -temp_im * im));
  y[hszCostab].im =
      0.5 * ((re + temp_im * y_re_tmp) + (temp2_re + -temp_im * b_re));
}

void FFTImplementationCallback::dobluesteinfft(
    const array<double, 1U> &x, int n2blue, int nfft,
    const array<double, 2U> &costab, const array<double, 2U> &sintab,
    const array<double, 2U> &sintabinv, array<creal_T, 1U> &y)
{
  array<creal_T, 1U> fv;
  array<creal_T, 1U> fy;
  array<creal_T, 1U> wwc;
  double nt_im;
  double nt_re;
  int ju;
  int minNrowsNx;
  int nInt2;
  int nInt2m1;
  int rt;
  unsigned int u;
  u = static_cast<unsigned int>(nfft) & 1U;
  if (u == 0U) {
    ju = static_cast<int>(static_cast<unsigned int>(nfft) >> 1);
    nInt2m1 = (ju + ju) - 1;
    wwc.set_size(nInt2m1);
    rt = 0;
    wwc[ju - 1].re = 1.0;
    wwc[ju - 1].im = 0.0;
    nInt2 = ju << 1;
    for (int k{0}; k <= ju - 2; k++) {
      minNrowsNx = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= minNrowsNx) {
        rt += minNrowsNx - nInt2;
      } else {
        rt += minNrowsNx;
      }
      nt_im = -3.1415926535897931 * static_cast<double>(rt) /
              static_cast<double>(ju);
      nt_re = std::cos(nt_im);
      nt_im = std::sin(nt_im);
      minNrowsNx = (ju - k) - 2;
      wwc[minNrowsNx].re = nt_re;
      wwc[minNrowsNx].im = -nt_im;
    }
    minNrowsNx = nInt2m1 - 1;
    for (int k{minNrowsNx}; k >= ju; k--) {
      wwc[k] = wwc[(nInt2m1 - k) - 1];
    }
  } else {
    ju = (nfft + nfft) - 1;
    wwc.set_size(ju);
    rt = 0;
    wwc[nfft - 1].re = 1.0;
    wwc[nfft - 1].im = 0.0;
    nInt2 = nfft << 1;
    for (int k{0}; k <= nfft - 2; k++) {
      minNrowsNx = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= minNrowsNx) {
        rt += minNrowsNx - nInt2;
      } else {
        rt += minNrowsNx;
      }
      nt_im = -3.1415926535897931 * static_cast<double>(rt) /
              static_cast<double>(nfft);
      nt_re = std::cos(nt_im);
      nt_im = std::sin(nt_im);
      minNrowsNx = (nfft - k) - 2;
      wwc[minNrowsNx].re = nt_re;
      wwc[minNrowsNx].im = -nt_im;
    }
    minNrowsNx = ju - 1;
    for (int k{minNrowsNx}; k >= nfft; k--) {
      wwc[k] = wwc[(ju - k) - 1];
    }
  }
  y.set_size(nfft);
  if (nfft > x.size(0)) {
    y.set_size(nfft);
    for (int k{0}; k < nfft; k++) {
      y[k].re = 0.0;
      y[k].im = 0.0;
    }
  }
  if ((n2blue != 1) && (u == 0U)) {
    FFTImplementationCallback::doHalfLengthBluestein(
        x, y, x.size(0), nfft, n2blue, wwc, costab, sintab, costab, sintabinv);
  } else {
    double im;
    double re;
    double twid_re;
    int b_k;
    int iDelta2;
    int iheight;
    int nRowsD2;
    minNrowsNx = x.size(0);
    if (nfft <= minNrowsNx) {
      minNrowsNx = nfft;
    }
    for (int k{0}; k < minNrowsNx; k++) {
      rt = (nfft + k) - 1;
      y[k].re = wwc[rt].re * x[k];
      y[k].im = wwc[rt].im * -x[k];
    }
    minNrowsNx++;
    for (int k{minNrowsNx}; k <= nfft; k++) {
      y[k - 1].re = 0.0;
      y[k - 1].im = 0.0;
    }
    fy.set_size(n2blue);
    if (n2blue > y.size(0)) {
      fy.set_size(n2blue);
      std::memset(&fy[0], 0,
                  static_cast<unsigned int>(n2blue) * sizeof(creal_T));
    }
    rt = y.size(0);
    if (rt > n2blue) {
      rt = n2blue;
    }
    nInt2m1 = n2blue - 2;
    nRowsD2 = static_cast<int>(static_cast<unsigned int>(n2blue) >> 1);
    b_k = static_cast<int>(static_cast<unsigned int>(nRowsD2) >> 1);
    nInt2 = 0;
    ju = 0;
    for (int k{0}; k <= rt - 2; k++) {
      bool tst;
      fy[nInt2] = y[k];
      minNrowsNx = n2blue;
      tst = true;
      while (tst) {
        minNrowsNx >>= 1;
        ju ^= minNrowsNx;
        tst = ((ju & minNrowsNx) == 0);
      }
      nInt2 = ju;
    }
    if (rt - 2 < 0) {
      minNrowsNx = 0;
    } else {
      minNrowsNx = rt - 1;
    }
    fy[nInt2] = y[minNrowsNx];
    if (n2blue > 1) {
      for (int k{0}; k <= nInt2m1; k += 2) {
        nt_re = fy[k + 1].re;
        nt_im = fy[k + 1].im;
        re = fy[k].re;
        im = fy[k].im;
        fy[k + 1].re = re - nt_re;
        fy[k + 1].im = fy[k].im - fy[k + 1].im;
        re += nt_re;
        im += nt_im;
        fy[k].re = re;
        fy[k].im = im;
      }
    }
    nInt2m1 = 2;
    iDelta2 = 4;
    iheight = ((b_k - 1) << 2) + 1;
    while (b_k > 0) {
      int i;
      for (i = 0; i < iheight; i += iDelta2) {
        minNrowsNx = i + nInt2m1;
        nt_re = fy[minNrowsNx].re;
        nt_im = fy[minNrowsNx].im;
        fy[minNrowsNx].re = fy[i].re - nt_re;
        fy[minNrowsNx].im = fy[i].im - nt_im;
        fy[i].re = fy[i].re + nt_re;
        fy[i].im = fy[i].im + nt_im;
      }
      minNrowsNx = 1;
      for (rt = b_k; rt < nRowsD2; rt += b_k) {
        double twid_im;
        twid_re = costab[rt];
        twid_im = sintab[rt];
        i = minNrowsNx;
        nInt2 = minNrowsNx + iheight;
        while (i < nInt2) {
          ju = i + nInt2m1;
          re = fy[ju].im;
          im = fy[ju].re;
          nt_re = twid_re * im - twid_im * re;
          nt_im = twid_re * re + twid_im * im;
          fy[ju].re = fy[i].re - nt_re;
          fy[ju].im = fy[i].im - nt_im;
          fy[i].re = fy[i].re + nt_re;
          fy[i].im = fy[i].im + nt_im;
          i += iDelta2;
        }
        minNrowsNx++;
      }
      b_k = static_cast<int>(static_cast<unsigned int>(b_k) >> 1);
      nInt2m1 = iDelta2;
      iDelta2 += iDelta2;
      iheight -= nInt2m1;
    }
    FFTImplementationCallback::r2br_r2dit_trig_impl(wwc, n2blue, costab, sintab,
                                                    fv);
    minNrowsNx = fy.size(0);
    for (int k{0}; k < minNrowsNx; k++) {
      re = fy[k].re;
      im = fv[k].im;
      nt_im = fy[k].im;
      twid_re = fv[k].re;
      fy[k].re = re * twid_re - nt_im * im;
      fy[k].im = re * im + nt_im * twid_re;
    }
    FFTImplementationCallback::r2br_r2dit_trig_impl(fy, n2blue, costab,
                                                    sintabinv, fv);
    if (fv.size(0) > 1) {
      re = 1.0 / static_cast<double>(fv.size(0));
      minNrowsNx = fv.size(0);
      for (int k{0}; k < minNrowsNx; k++) {
        fv[k].re = re * fv[k].re;
        fv[k].im = re * fv[k].im;
      }
    }
    minNrowsNx = wwc.size(0);
    for (int k{nfft}; k <= minNrowsNx; k++) {
      re = wwc[k - 1].re;
      im = fv[k - 1].im;
      nt_im = wwc[k - 1].im;
      twid_re = fv[k - 1].re;
      rt = k - nfft;
      y[rt].re = re * twid_re + nt_im * im;
      y[rt].im = re * im - nt_im * twid_re;
    }
  }
}

} // namespace fft
} // namespace internal
} // namespace coder

// End of code generation (FFTImplementationCallback.cpp)
