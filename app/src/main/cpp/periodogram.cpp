//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// periodogram.cpp
//
// Code generation for function 'periodogram'
//

// Include files
#include "periodogram.h"
#include "FFTImplementationCallback.h"
#include "rt_nonfinite.h"
#include "stickyStruct.h"
#include "coder_array.h"
#include <algorithm>
#include <cmath>
#include <cstring>

// Function Definitions
namespace coder {
void periodogram(const array<double, 1U> &x, array<double, 1U> &Px,
                 array<double, 1U> &w)
{
  internal::i_stickyStruct obj;
  array<creal_T, 1U> Xx;
  array<double, 2U> b_select;
  array<double, 2U> costab1q;
  array<double, 2U> sintab;
  array<double, 2U> sintabinv;
  array<double, 2U> wrappedData;
  array<double, 1U> b_xw;
  array<double, 1U> xw;
  double b_halfNPTS_tmp;
  double freq_res;
  double halfNPTS_tmp;
  double half_res;
  double isNPTSodd_tmp;
  int acoef;
  int csz_idx_0;
  int loop_ub;
  int nFullPasses;
  bool useRadix2;
  freq_res = std::frexp(static_cast<double>(x.size(0)), &csz_idx_0);
  if (freq_res == 0.5) {
    csz_idx_0--;
  }
  if (csz_idx_0 == 0) {
    freq_res = 1.0;
  } else if (csz_idx_0 == 1) {
    freq_res = 2.0;
  } else if (csz_idx_0 == 2) {
    freq_res = 4.0;
  } else {
    freq_res = std::pow(2.0, static_cast<double>(csz_idx_0));
  }
  obj.next.next.next.next.next.next.next.next.value =
      std::fmax(256.0, freq_res);
  if (x.size(0) == 1) {
    csz_idx_0 = x.size(0);
  } else if (x.size(0) == 1) {
    csz_idx_0 = 1;
  } else {
    csz_idx_0 = x.size(0);
  }
  xw.set_size(csz_idx_0);
  if (csz_idx_0 != 0) {
    acoef = (x.size(0) != 1);
    for (int k{0}; k < csz_idx_0; k++) {
      xw[k] = x[acoef * k];
    }
  }
  if (xw.size(0) > obj.next.next.next.next.next.next.next.next.value) {
    loop_ub =
        static_cast<int>(obj.next.next.next.next.next.next.next.next.value);
    b_xw.set_size(
        static_cast<int>(obj.next.next.next.next.next.next.next.next.value));
    std::memset(&b_xw[0], 0,
                static_cast<unsigned int>(loop_ub) * sizeof(double));
    wrappedData.set_size(
        static_cast<int>(obj.next.next.next.next.next.next.next.next.value), 1);
    std::memset(&wrappedData[0], 0,
                static_cast<unsigned int>(loop_ub) * sizeof(double));
    nFullPasses = static_cast<int>(
        static_cast<unsigned int>(xw.size(0)) /
        static_cast<unsigned int>(static_cast<int>(
            obj.next.next.next.next.next.next.next.next.value)));
    csz_idx_0 =
        nFullPasses *
        static_cast<int>(obj.next.next.next.next.next.next.next.next.value);
    acoef = (xw.size(0) - csz_idx_0) - 1;
    for (int k{0}; k <= acoef; k++) {
      wrappedData[k] = xw[csz_idx_0 + k];
    }
    csz_idx_0 = acoef + 2;
    if (csz_idx_0 <= loop_ub) {
      std::memset(&wrappedData[csz_idx_0 + -1], 0,
                  static_cast<unsigned int>((loop_ub - csz_idx_0) + 1) *
                      sizeof(double));
    }
    for (int k{0}; k < nFullPasses; k++) {
      csz_idx_0 = k * static_cast<int>(
                          obj.next.next.next.next.next.next.next.next.value);
      for (int b_k{0}; b_k < loop_ub; b_k++) {
        wrappedData[b_k] = wrappedData[b_k] + xw[csz_idx_0 + b_k];
      }
    }
    std::copy(&wrappedData[0], &wrappedData[loop_ub], &b_xw[0]);
  } else {
    b_xw.set_size(csz_idx_0);
    if (csz_idx_0 - 1 >= 0) {
      std::copy(&xw[0], &xw[csz_idx_0], &b_xw[0]);
    }
  }
  if (b_xw.size(0) == 0) {
    csz_idx_0 =
        static_cast<int>(obj.next.next.next.next.next.next.next.next.value);
    Xx.set_size(
        static_cast<int>(obj.next.next.next.next.next.next.next.next.value));
    if (csz_idx_0 - 1 >= 0) {
      std::memset(&Xx[0], 0,
                  static_cast<unsigned int>(csz_idx_0) * sizeof(creal_T));
    }
  } else {
    int n2;
    useRadix2 =
        ((static_cast<int>(obj.next.next.next.next.next.next.next.next.value) &
          (static_cast<int>(obj.next.next.next.next.next.next.next.next.value) -
           1)) == 0);
    loop_ub = 1;
    if (useRadix2) {
      csz_idx_0 =
          static_cast<int>(obj.next.next.next.next.next.next.next.next.value);
    } else {
      bool exitg1;
      csz_idx_0 =
          (static_cast<int>(obj.next.next.next.next.next.next.next.next.value) +
           static_cast<int>(
               obj.next.next.next.next.next.next.next.next.value)) -
          1;
      acoef = 31;
      n2 = 0;
      exitg1 = false;
      while ((!exitg1) && (acoef - n2 > 1)) {
        loop_ub = (n2 + acoef) >> 1;
        nFullPasses = 1 << loop_ub;
        if (nFullPasses == csz_idx_0) {
          acoef = loop_ub;
          exitg1 = true;
        } else if (nFullPasses > csz_idx_0) {
          acoef = loop_ub;
        } else {
          n2 = loop_ub;
        }
      }
      loop_ub = 1 << acoef;
      csz_idx_0 = loop_ub;
    }
    freq_res = 6.2831853071795862 / static_cast<double>(csz_idx_0);
    acoef = static_cast<int>((static_cast<unsigned int>(csz_idx_0) >> 1) >> 1);
    costab1q.set_size(1, acoef + 1);
    costab1q[0] = 1.0;
    csz_idx_0 = static_cast<int>(static_cast<unsigned int>(acoef) >> 1) - 1;
    for (int k{0}; k <= csz_idx_0; k++) {
      costab1q[k + 1] = std::cos(freq_res * (static_cast<double>(k) + 1.0));
    }
    csz_idx_0 += 2;
    for (int k{csz_idx_0}; k < acoef; k++) {
      costab1q[k] = std::sin(freq_res * static_cast<double>(acoef - k));
    }
    costab1q[acoef] = 0.0;
    if (!useRadix2) {
      acoef = costab1q.size(1) - 1;
      n2 = (costab1q.size(1) - 1) << 1;
      b_select.set_size(1, n2 + 1);
      sintab.set_size(1, n2 + 1);
      b_select[0] = 1.0;
      sintab[0] = 0.0;
      sintabinv.set_size(1, n2 + 1);
      for (int k{0}; k < acoef; k++) {
        sintabinv[k + 1] = costab1q[(acoef - k) - 1];
      }
      csz_idx_0 = costab1q.size(1);
      for (int k{csz_idx_0}; k <= n2; k++) {
        sintabinv[k] = costab1q[k - acoef];
      }
      for (int k{0}; k < acoef; k++) {
        b_select[k + 1] = costab1q[k + 1];
        sintab[k + 1] = -costab1q[(acoef - k) - 1];
      }
      for (int k{csz_idx_0}; k <= n2; k++) {
        b_select[k] = -costab1q[n2 - k];
        sintab[k] = -costab1q[k - acoef];
      }
      internal::fft::FFTImplementationCallback::dobluesteinfft(
          b_xw, loop_ub,
          static_cast<int>(obj.next.next.next.next.next.next.next.next.value),
          b_select, sintab, sintabinv, Xx);
    } else {
      acoef = costab1q.size(1) - 1;
      n2 = (costab1q.size(1) - 1) << 1;
      b_select.set_size(1, n2 + 1);
      sintab.set_size(1, n2 + 1);
      b_select[0] = 1.0;
      sintab[0] = 0.0;
      for (int k{0}; k < acoef; k++) {
        b_select[k + 1] = costab1q[k + 1];
        sintab[k + 1] = -costab1q[(acoef - k) - 1];
      }
      csz_idx_0 = costab1q.size(1);
      for (int k{csz_idx_0}; k <= n2; k++) {
        b_select[k] = -costab1q[n2 - k];
        sintab[k] = -costab1q[k - acoef];
      }
      csz_idx_0 =
          static_cast<int>(obj.next.next.next.next.next.next.next.next.value);
      Xx.set_size(
          static_cast<int>(obj.next.next.next.next.next.next.next.next.value));
      if (static_cast<int>(obj.next.next.next.next.next.next.next.next.value) >
          b_xw.size(0)) {
        Xx.set_size(static_cast<int>(
            obj.next.next.next.next.next.next.next.next.value));
        if (csz_idx_0 - 1 >= 0) {
          std::memset(&Xx[0], 0,
                      static_cast<unsigned int>(csz_idx_0) * sizeof(creal_T));
        }
      }
      internal::fft::FFTImplementationCallback::doHalfLengthRadix2(
          b_xw, Xx,
          static_cast<int>(obj.next.next.next.next.next.next.next.next.value),
          b_select, sintab);
    }
  }
  freq_res = 50.0 / obj.next.next.next.next.next.next.next.next.value;
  costab1q.set_size(
      1, static_cast<int>(obj.next.next.next.next.next.next.next.next.value -
                          1.0) +
             1);
  csz_idx_0 =
      static_cast<int>(obj.next.next.next.next.next.next.next.next.value - 1.0);
  for (int k{0}; k <= csz_idx_0; k++) {
    costab1q[k] = k;
  }
  costab1q.set_size(1, costab1q.size(1));
  csz_idx_0 = costab1q.size(1) - 1;
  for (int k{0}; k <= csz_idx_0; k++) {
    costab1q[k] = freq_res * costab1q[k];
  }
  half_res = freq_res / 2.0;
  if (std::isinf(obj.next.next.next.next.next.next.next.next.value)) {
    isNPTSodd_tmp = rtNaN;
  } else {
    isNPTSodd_tmp =
        std::fmod(obj.next.next.next.next.next.next.next.next.value, 2.0);
  }
  useRadix2 = (isNPTSodd_tmp != 0.0);
  halfNPTS_tmp = obj.next.next.next.next.next.next.next.next.value / 2.0 + 1.0;
  b_halfNPTS_tmp =
      (obj.next.next.next.next.next.next.next.next.value + 1.0) / 2.0;
  if (useRadix2) {
    costab1q[static_cast<int>(b_halfNPTS_tmp) - 1] = 25.0 - half_res;
    costab1q[static_cast<int>(static_cast<unsigned int>(b_halfNPTS_tmp))] =
        half_res + 25.0;
  } else {
    costab1q[static_cast<int>(halfNPTS_tmp) - 1] = 25.0;
  }
  costab1q[static_cast<int>(obj.next.next.next.next.next.next.next.next.value) -
           1] = 50.0 - freq_res;
  csz_idx_0 = x.size(0);
  acoef = Xx.size(0);
  xw.set_size(Xx.size(0));
  for (int k{0}; k < acoef; k++) {
    double Xx_re;
    freq_res = Xx[k].re;
    half_res = Xx[k].im;
    Xx_re = freq_res * freq_res - half_res * -half_res;
    if (freq_res * -half_res + half_res * freq_res == 0.0) {
      freq_res = Xx_re / static_cast<double>(csz_idx_0);
    } else if (Xx_re == 0.0) {
      freq_res = 0.0;
    } else {
      freq_res = Xx_re / static_cast<double>(csz_idx_0);
    }
    xw[k] = freq_res;
  }
  if (isNPTSodd_tmp == 1.0) {
    b_select.set_size(1, static_cast<int>(b_halfNPTS_tmp - 1.0) + 1);
    csz_idx_0 = static_cast<int>(b_halfNPTS_tmp - 1.0);
    for (int k{0}; k <= csz_idx_0; k++) {
      b_select[k] = static_cast<double>(k) + 1.0;
    }
    Px.set_size(static_cast<int>(b_halfNPTS_tmp - 1.0) + 1);
    for (int k{0}; k <= csz_idx_0; k++) {
      Px[k] = xw[k];
    }
    b_xw.set_size(static_cast<int>(b_halfNPTS_tmp - 1.0));
    for (int k{0}; k < csz_idx_0; k++) {
      b_xw[k] = 2.0 * Px[k + 1];
    }
    for (int k{0}; k < csz_idx_0; k++) {
      Px[k + 1] = b_xw[k];
    }
  } else {
    b_select.set_size(1, static_cast<int>(halfNPTS_tmp - 1.0) + 1);
    acoef = static_cast<int>(halfNPTS_tmp - 1.0);
    for (int k{0}; k <= acoef; k++) {
      b_select[k] = static_cast<double>(k) + 1.0;
    }
    Px.set_size(static_cast<int>(halfNPTS_tmp - 1.0) + 1);
    for (int k{0}; k <= acoef; k++) {
      Px[k] = xw[k];
    }
    csz_idx_0 = static_cast<int>(halfNPTS_tmp - 1.0) - 1;
    b_xw.set_size(static_cast<int>(halfNPTS_tmp - 1.0) - 1);
    for (int k{0}; k <= acoef - 2; k++) {
      b_xw[k] = 2.0 * Px[k + 1];
    }
    for (int k{0}; k < csz_idx_0; k++) {
      Px[k + 1] = b_xw[k];
    }
  }
  csz_idx_0 = b_select.size(1);
  w.set_size(b_select.size(1));
  for (int k{0}; k < csz_idx_0; k++) {
    w[k] = costab1q[static_cast<int>(b_select[k]) - 1];
  }
  for (int k{0}; k < csz_idx_0; k++) {
    Px[k] = Px[k] / 50.0;
  }
}

} // namespace coder

// End of code generation (periodogram.cpp)
