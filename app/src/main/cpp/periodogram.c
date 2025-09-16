/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * periodogram.c
 *
 * Code generation for function 'periodogram'
 *
 */

/* Include files */
#include "periodogram.h"
#include "FFTImplementationCallback.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "omp.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>

/* Function Definitions */
void periodogram(const double x_data[], int x_size, double varargin_3,
                 emxArray_real_T *Px, emxArray_real_T *w)
{
  __m128d r;
  emxArray_creal_T *Xx;
  emxArray_real_T *Sxx;
  emxArray_real_T *b_select;
  emxArray_real_T *costab1q;
  emxArray_real_T *sintab;
  emxArray_real_T *sintabinv;
  creal_T x;
  creal_T *Xx_data;
  double xw_data[150];
  double dv[2];
  double Fs;
  double Nyq;
  double Xx_im_tmp;
  double Xx_re;
  double b_Xx_re;
  double b_halfNPTS_tmp;
  double f;
  double freq_res;
  double halfNPTS_tmp;
  double half_res;
  double isNPTSodd_tmp;
  double obj_next_next_value;
  double *Px_data;
  double *costab1q_data;
  double *select_data;
  double *w_data;
  int acoef;
  int b_k;
  int b_n2;
  int c_k;
  int csz_idx_0;
  int d_k;
  int e_k;
  int f_k;
  int g_k;
  int h_k;
  int i;
  int k;
  int loop_ub;
  int n;
  int n2;
  boolean_T Fs_tmp;
  boolean_T useRadix2;
  f = frexp(x_size, &acoef);
  if (f == 0.5) {
    acoef--;
  }
  if (acoef == 0) {
    f = 1.0;
  } else if (acoef == 1) {
    f = 2.0;
  } else if (acoef == 2) {
    f = 4.0;
  } else {
    f = pow(2.0, acoef);
  }
  obj_next_next_value = fmax(256.0, f);
  Fs_tmp = rtIsNaN(varargin_3);
  if (Fs_tmp) {
    Fs = 6.2831853071795862;
  } else {
    Fs = varargin_3;
  }
  if (x_size == 1) {
    csz_idx_0 = x_size;
  } else if (x_size == 1) {
    csz_idx_0 = 1;
  } else {
    csz_idx_0 = x_size;
  }
  if (csz_idx_0 != 0) {
    acoef = (x_size != 1);
    for (k = 0; k < csz_idx_0; k++) {
      xw_data[k] = x_data[acoef * k];
    }
  }
  emxInit_real_T(&b_select, 2);
  emxInit_creal_T(&Xx);
  emxInit_real_T(&costab1q, 2);
  if (csz_idx_0 == 0) {
    n2 = (int)obj_next_next_value;
    acoef = Xx->size[0];
    Xx->size[0] = (int)obj_next_next_value;
    emxEnsureCapacity_creal_T(Xx, acoef);
    Xx_data = Xx->data;
    for (k = 0; k < n2; k++) {
      Xx_data[k].re = 0.0;
      Xx_data[k].im = 0.0;
    }
  } else {
    int N2blue;
    useRadix2 =
        (((int)obj_next_next_value & ((int)obj_next_next_value - 1)) == 0);
    N2blue = 1;
    if (useRadix2) {
      acoef = (int)obj_next_next_value;
    } else {
      boolean_T exitg1;
      acoef = ((int)obj_next_next_value + (int)obj_next_next_value) - 1;
      n2 = 31;
      b_n2 = 0;
      exitg1 = false;
      while ((!exitg1) && (n2 - b_n2 > 1)) {
        loop_ub = (b_n2 + n2) >> 1;
        n = 1 << loop_ub;
        if (n == acoef) {
          n2 = loop_ub;
          exitg1 = true;
        } else if (n > acoef) {
          n2 = loop_ub;
        } else {
          b_n2 = loop_ub;
        }
      }
      N2blue = 1 << n2;
      acoef = N2blue;
    }
    f = 6.2831853071795862 / (double)acoef;
    b_n2 = (int)(((unsigned int)acoef >> 1) >> 1);
    acoef = costab1q->size[0] * costab1q->size[1];
    costab1q->size[0] = 1;
    costab1q->size[1] = b_n2 + 1;
    emxEnsureCapacity_real_T(costab1q, acoef);
    costab1q_data = costab1q->data;
    costab1q_data[0] = 1.0;
    acoef = (int)((unsigned int)b_n2 >> 1) - 1;
    if (acoef + 1 < 1200) {
      for (b_k = 0; b_k <= acoef; b_k++) {
        costab1q_data[b_k + 1] = cos(f * ((double)b_k + 1.0));
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

      for (b_k = 0; b_k <= acoef; b_k++) {
        costab1q_data[b_k + 1] = cos(f * ((double)b_k + 1.0));
      }
    }
    n2 = acoef + 2;
    if ((b_n2 - acoef) - 2 < 1200) {
      for (c_k = n2; c_k < b_n2; c_k++) {
        costab1q_data[c_k] = sin(f * (double)(b_n2 - c_k));
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

      for (c_k = n2; c_k < b_n2; c_k++) {
        costab1q_data[c_k] = sin(f * (double)(b_n2 - c_k));
      }
    }
    costab1q_data[b_n2] = 0.0;
    emxInit_real_T(&sintab, 2);
    emxInit_real_T(&sintabinv, 2);
    if (!useRadix2) {
      n = costab1q->size[1] - 1;
      n2 = (costab1q->size[1] - 1) << 1;
      acoef = b_select->size[0] * b_select->size[1];
      b_select->size[0] = 1;
      b_select->size[1] = n2 + 1;
      emxEnsureCapacity_real_T(b_select, acoef);
      select_data = b_select->data;
      acoef = sintab->size[0] * sintab->size[1];
      sintab->size[0] = 1;
      sintab->size[1] = n2 + 1;
      emxEnsureCapacity_real_T(sintab, acoef);
      w_data = sintab->data;
      select_data[0] = 1.0;
      w_data[0] = 0.0;
      acoef = sintabinv->size[0] * sintabinv->size[1];
      sintabinv->size[0] = 1;
      sintabinv->size[1] = n2 + 1;
      emxEnsureCapacity_real_T(sintabinv, acoef);
      Px_data = sintabinv->data;
      acoef = (costab1q->size[1] - 1 < 1200);
      if (acoef) {
        for (f_k = 0; f_k < n; f_k++) {
          Px_data[f_k + 1] = costab1q_data[(n - f_k) - 1];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (f_k = 0; f_k < n; f_k++) {
          Px_data[f_k + 1] = costab1q_data[(n - f_k) - 1];
        }
      }
      loop_ub = costab1q->size[1];
      for (k = loop_ub; k <= n2; k++) {
        Px_data[k] = costab1q_data[k - n];
      }
      if (acoef) {
        for (g_k = 0; g_k < n; g_k++) {
          select_data[g_k + 1] = costab1q_data[g_k + 1];
          w_data[g_k + 1] = -costab1q_data[(n - g_k) - 1];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (g_k = 0; g_k < n; g_k++) {
          select_data[g_k + 1] = costab1q_data[g_k + 1];
          w_data[g_k + 1] = -costab1q_data[(n - g_k) - 1];
        }
      }
      if ((n2 - costab1q->size[1]) + 1 < 1200) {
        for (h_k = loop_ub; h_k <= n2; h_k++) {
          select_data[h_k] = -costab1q_data[n2 - h_k];
          w_data[h_k] = -costab1q_data[h_k - n];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (h_k = loop_ub; h_k <= n2; h_k++) {
          select_data[h_k] = -costab1q_data[n2 - h_k];
          w_data[h_k] = -costab1q_data[h_k - n];
        }
      }
      c_FFTImplementationCallback_dob(xw_data, csz_idx_0, N2blue,
                                      (int)obj_next_next_value, b_select,
                                      sintab, sintabinv, Xx);
      Xx_data = Xx->data;
    } else {
      n2 = costab1q->size[1] - 1;
      b_n2 = (costab1q->size[1] - 1) << 1;
      acoef = b_select->size[0] * b_select->size[1];
      b_select->size[0] = 1;
      b_select->size[1] = b_n2 + 1;
      emxEnsureCapacity_real_T(b_select, acoef);
      select_data = b_select->data;
      acoef = sintab->size[0] * sintab->size[1];
      sintab->size[0] = 1;
      sintab->size[1] = b_n2 + 1;
      emxEnsureCapacity_real_T(sintab, acoef);
      w_data = sintab->data;
      select_data[0] = 1.0;
      w_data[0] = 0.0;
      if (costab1q->size[1] - 1 < 1200) {
        for (d_k = 0; d_k < n2; d_k++) {
          select_data[d_k + 1] = costab1q_data[d_k + 1];
          w_data[d_k + 1] = -costab1q_data[(n2 - d_k) - 1];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (d_k = 0; d_k < n2; d_k++) {
          select_data[d_k + 1] = costab1q_data[d_k + 1];
          w_data[d_k + 1] = -costab1q_data[(n2 - d_k) - 1];
        }
      }
      acoef = costab1q->size[1];
      if ((b_n2 - costab1q->size[1]) + 1 < 1200) {
        for (e_k = acoef; e_k <= b_n2; e_k++) {
          select_data[e_k] = -costab1q_data[b_n2 - e_k];
          w_data[e_k] = -costab1q_data[e_k - n2];
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (e_k = acoef; e_k <= b_n2; e_k++) {
          select_data[e_k] = -costab1q_data[b_n2 - e_k];
          w_data[e_k] = -costab1q_data[e_k - n2];
        }
      }
      n2 = (int)obj_next_next_value;
      acoef = Xx->size[0];
      Xx->size[0] = (int)obj_next_next_value;
      emxEnsureCapacity_creal_T(Xx, acoef);
      if ((int)obj_next_next_value > csz_idx_0) {
        acoef = Xx->size[0];
        Xx->size[0] = (int)obj_next_next_value;
        emxEnsureCapacity_creal_T(Xx, acoef);
        Xx_data = Xx->data;
        for (k = 0; k < n2; k++) {
          Xx_data[k].re = 0.0;
          Xx_data[k].im = 0.0;
        }
      }
      c_FFTImplementationCallback_doH(
          xw_data, csz_idx_0, Xx, (int)obj_next_next_value, b_select, sintab);
      Xx_data = Xx->data;
    }
    emxFree_real_T(&sintabinv);
    emxFree_real_T(&sintab);
  }
  freq_res = Fs / obj_next_next_value;
  acoef = costab1q->size[0] * costab1q->size[1];
  costab1q->size[0] = 1;
  costab1q->size[1] = (int)(obj_next_next_value - 1.0) + 1;
  emxEnsureCapacity_real_T(costab1q, acoef);
  costab1q_data = costab1q->data;
  acoef = (int)(obj_next_next_value - 1.0);
  for (k = 0; k <= acoef; k++) {
    costab1q_data[k] = k;
  }
  acoef = costab1q->size[0] * costab1q->size[1];
  costab1q->size[0] = 1;
  emxEnsureCapacity_real_T(costab1q, acoef);
  costab1q_data = costab1q->data;
  acoef = costab1q->size[1] - 1;
  n2 = (costab1q->size[1] / 2) << 1;
  b_n2 = n2 - 2;
  for (k = 0; k <= b_n2; k += 2) {
    r = _mm_loadu_pd(&costab1q_data[k]);
    _mm_storeu_pd(&costab1q_data[k], _mm_mul_pd(_mm_set1_pd(freq_res), r));
  }
  for (k = n2; k <= acoef; k++) {
    costab1q_data[k] *= freq_res;
  }
  Nyq = Fs / 2.0;
  half_res = freq_res / 2.0;
  isNPTSodd_tmp = fmod(obj_next_next_value, 2.0);
  if (rtIsInf(obj_next_next_value)) {
    f = rtNaN;
  } else {
    f = isNPTSodd_tmp;
  }
  useRadix2 = (f != 0.0);
  halfNPTS_tmp = obj_next_next_value / 2.0 + 1.0;
  b_halfNPTS_tmp = (obj_next_next_value + 1.0) / 2.0;
  if (useRadix2) {
    costab1q_data[(int)b_halfNPTS_tmp - 1] = Nyq - half_res;
    costab1q_data[(int)(unsigned int)b_halfNPTS_tmp] = Nyq + half_res;
  } else {
    costab1q_data[(int)halfNPTS_tmp - 1] = Nyq;
  }
  costab1q_data[(int)obj_next_next_value - 1] = Fs - freq_res;
  x.re = x_size;
  emxInit_real_T(&Sxx, 1);
  n2 = Xx->size[0];
  acoef = Sxx->size[0];
  Sxx->size[0] = Xx->size[0];
  emxEnsureCapacity_real_T(Sxx, acoef);
  w_data = Sxx->data;
  acoef = Xx->size[0];
  if (Xx->size[0] < 1200) {
    for (i = 0; i < n2; i++) {
      f = Xx_data[i].re;
      freq_res = Xx_data[i].im;
      Nyq = f * f - freq_res * -freq_res;
      if (f * -freq_res + freq_res * f == 0.0) {
        f = Nyq / (double)x_size;
      } else if (Nyq == 0.0) {
        f = 0.0;
      } else {
        f = Nyq / (double)x_size;
      }
      w_data[i] = f;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        Xx_re, Xx_im_tmp, b_Xx_re)

    for (i = 0; i < acoef; i++) {
      Xx_re = Xx_data[i].re;
      Xx_im_tmp = Xx_data[i].im;
      b_Xx_re = Xx_re * Xx_re - Xx_im_tmp * -Xx_im_tmp;
      if (Xx_re * -Xx_im_tmp + Xx_im_tmp * Xx_re == 0.0) {
        Xx_re = b_Xx_re / x.re;
      } else if (b_Xx_re == 0.0) {
        Xx_re = 0.0;
      } else {
        Xx_re = b_Xx_re / x.re;
      }
      w_data[i] = Xx_re;
    }
  }
  emxFree_creal_T(&Xx);
  if ((int)isNPTSodd_tmp == 1) {
    acoef = b_select->size[0] * b_select->size[1];
    b_select->size[0] = 1;
    b_select->size[1] = (int)(b_halfNPTS_tmp - 1.0) + 1;
    emxEnsureCapacity_real_T(b_select, acoef);
    select_data = b_select->data;
    loop_ub = (int)(b_halfNPTS_tmp - 1.0);
    acoef = (((int)(b_halfNPTS_tmp - 1.0) + 1) / 2) << 1;
    n2 = acoef - 2;
    for (k = 0; k <= n2; k += 2) {
      dv[0] = k;
      dv[1] = k + 1;
      r = _mm_loadu_pd(&dv[0]);
      _mm_storeu_pd(&select_data[k], _mm_add_pd(_mm_set1_pd(1.0), r));
    }
    for (k = acoef; k <= loop_ub; k++) {
      select_data[k] = (double)k + 1.0;
    }
    acoef = Px->size[0];
    Px->size[0] = (int)(b_halfNPTS_tmp - 1.0) + 1;
    emxEnsureCapacity_real_T(Px, acoef);
    Px_data = Px->data;
    for (k = 0; k <= loop_ub; k++) {
      Px_data[k] = w_data[k];
    }
    acoef = Sxx->size[0];
    Sxx->size[0] = (int)(b_halfNPTS_tmp - 1.0);
    emxEnsureCapacity_real_T(Sxx, acoef);
    w_data = Sxx->data;
    acoef = ((int)(b_halfNPTS_tmp - 1.0) / 2) << 1;
    n2 = acoef - 2;
    for (k = 0; k <= n2; k += 2) {
      r = _mm_loadu_pd(&Px_data[k + 1]);
      _mm_storeu_pd(&w_data[k], _mm_mul_pd(_mm_set1_pd(2.0), r));
    }
    for (k = acoef; k < loop_ub; k++) {
      w_data[k] = 2.0 * Px_data[k + 1];
    }
    acoef = Sxx->size[0];
    for (k = 0; k < acoef; k++) {
      Px_data[k + 1] = w_data[k];
    }
  } else {
    acoef = b_select->size[0] * b_select->size[1];
    b_select->size[0] = 1;
    b_select->size[1] = (int)(halfNPTS_tmp - 1.0) + 1;
    emxEnsureCapacity_real_T(b_select, acoef);
    select_data = b_select->data;
    n = (int)(halfNPTS_tmp - 1.0);
    acoef = (((int)(halfNPTS_tmp - 1.0) + 1) / 2) << 1;
    n2 = acoef - 2;
    for (k = 0; k <= n2; k += 2) {
      dv[0] = k;
      dv[1] = k + 1;
      r = _mm_loadu_pd(&dv[0]);
      _mm_storeu_pd(&select_data[k], _mm_add_pd(_mm_set1_pd(1.0), r));
    }
    for (k = acoef; k <= n; k++) {
      select_data[k] = (double)k + 1.0;
    }
    acoef = Px->size[0];
    Px->size[0] = (int)(halfNPTS_tmp - 1.0) + 1;
    emxEnsureCapacity_real_T(Px, acoef);
    Px_data = Px->data;
    for (k = 0; k <= n; k++) {
      Px_data[k] = w_data[k];
    }
    acoef = Sxx->size[0];
    Sxx->size[0] = (int)(halfNPTS_tmp - 1.0) - 1;
    emxEnsureCapacity_real_T(Sxx, acoef);
    w_data = Sxx->data;
    acoef = (((int)(halfNPTS_tmp - 1.0) - 1) / 2) << 1;
    n2 = acoef - 2;
    for (k = 0; k <= n2; k += 2) {
      r = _mm_loadu_pd(&Px_data[k + 1]);
      _mm_storeu_pd(&w_data[k], _mm_mul_pd(_mm_set1_pd(2.0), r));
    }
    for (k = acoef; k <= n - 2; k++) {
      w_data[k] = 2.0 * Px_data[k + 1];
    }
    acoef = Sxx->size[0];
    for (k = 0; k < acoef; k++) {
      Px_data[k + 1] = w_data[k];
    }
  }
  emxFree_real_T(&Sxx);
  n2 = b_select->size[1];
  acoef = w->size[0];
  w->size[0] = b_select->size[1];
  emxEnsureCapacity_real_T(w, acoef);
  w_data = w->data;
  for (k = 0; k < n2; k++) {
    w_data[k] = costab1q_data[(int)select_data[k] - 1];
  }
  emxFree_real_T(&costab1q);
  emxFree_real_T(&b_select);
  if (!Fs_tmp) {
    acoef = Px->size[0];
    n2 = (Px->size[0] / 2) << 1;
    b_n2 = n2 - 2;
    for (k = 0; k <= b_n2; k += 2) {
      r = _mm_loadu_pd(&Px_data[k]);
      _mm_storeu_pd(&Px_data[k], _mm_div_pd(r, _mm_set1_pd(varargin_3)));
    }
    for (k = n2; k < acoef; k++) {
      Px_data[k] /= varargin_3;
    }
  } else {
    acoef = Px->size[0];
    n2 = (Px->size[0] / 2) << 1;
    b_n2 = n2 - 2;
    for (k = 0; k <= b_n2; k += 2) {
      r = _mm_loadu_pd(&Px_data[k]);
      _mm_storeu_pd(&Px_data[k],
                    _mm_div_pd(r, _mm_set1_pd(6.2831853071795862)));
    }
    for (k = n2; k < acoef; k++) {
      Px_data[k] /= 6.2831853071795862;
    }
  }
}

/* End of code generation (periodogram.c) */
