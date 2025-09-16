/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * FFTImplementationCallback.c
 *
 * Code generation for function 'FFTImplementationCallback'
 *
 */

/* Include files */
#include "FFTImplementationCallback.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "omp.h"
#include <math.h>

/* Function Declarations */
static void c_FFTImplementationCallback_r2b(const emxArray_creal_T *x,
                                            int unsigned_nRows,
                                            const emxArray_real_T *costab,
                                            const emxArray_real_T *sintab,
                                            emxArray_creal_T *y);

static void d_FFTImplementationCallback_doH(
    const double x_data[], int x_size, emxArray_creal_T *y, int nrowsx,
    int nRows, int nfft, const emxArray_creal_T *wwc,
    const emxArray_real_T *costab, const emxArray_real_T *sintab,
    const emxArray_real_T *costabinv, const emxArray_real_T *sintabinv);

/* Function Definitions */
static void c_FFTImplementationCallback_r2b(const emxArray_creal_T *x,
                                            int unsigned_nRows,
                                            const emxArray_real_T *costab,
                                            const emxArray_real_T *sintab,
                                            emxArray_creal_T *y)
{
  const creal_T *x_data;
  creal_T *y_data;
  const double *costab_data;
  const double *sintab_data;
  double im;
  double re;
  double temp_im;
  double temp_re;
  int b_i;
  int i;
  int iDelta;
  int iDelta2;
  int iheight;
  int istart;
  int iy;
  int j;
  int ju;
  int k;
  int nRowsD2;
  sintab_data = sintab->data;
  costab_data = costab->data;
  x_data = x->data;
  istart = y->size[0];
  y->size[0] = unsigned_nRows;
  emxEnsureCapacity_creal_T(y, istart);
  y_data = y->data;
  if (unsigned_nRows > x->size[0]) {
    istart = y->size[0];
    y->size[0] = unsigned_nRows;
    emxEnsureCapacity_creal_T(y, istart);
    y_data = y->data;
    for (i = 0; i < unsigned_nRows; i++) {
      y_data[i].re = 0.0;
      y_data[i].im = 0.0;
    }
  }
  j = x->size[0];
  if (j > unsigned_nRows) {
    j = unsigned_nRows;
  }
  iDelta = unsigned_nRows - 2;
  nRowsD2 = (int)((unsigned int)unsigned_nRows >> 1);
  k = (int)((unsigned int)nRowsD2 >> 1);
  iy = 0;
  ju = 0;
  for (i = 0; i <= j - 2; i++) {
    boolean_T tst;
    y_data[iy] = x_data[i];
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
  y_data[iy] = x_data[istart];
  if (unsigned_nRows > 1) {
    for (i = 0; i <= iDelta; i += 2) {
      temp_re = y_data[i + 1].re;
      temp_im = y_data[i + 1].im;
      re = y_data[i].re;
      im = y_data[i].im;
      y_data[i + 1].re = re - temp_re;
      y_data[i + 1].im = y_data[i].im - y_data[i + 1].im;
      re += temp_re;
      im += temp_im;
      y_data[i].re = re;
      y_data[i].im = im;
    }
  }
  iDelta = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    for (b_i = 0; b_i < iheight; b_i += iDelta2) {
      istart = b_i + iDelta;
      temp_re = y_data[istart].re;
      temp_im = y_data[istart].im;
      y_data[istart].re = y_data[b_i].re - temp_re;
      y_data[istart].im = y_data[b_i].im - temp_im;
      y_data[b_i].re += temp_re;
      y_data[b_i].im += temp_im;
    }
    istart = 1;
    for (j = k; j < nRowsD2; j += k) {
      double twid_im;
      double twid_re;
      twid_re = costab_data[j];
      twid_im = sintab_data[j];
      b_i = istart;
      iy = istart + iheight;
      while (b_i < iy) {
        ju = b_i + iDelta;
        re = y_data[ju].im;
        im = y_data[ju].re;
        temp_re = twid_re * im - twid_im * re;
        temp_im = twid_re * re + twid_im * im;
        y_data[ju].re = y_data[b_i].re - temp_re;
        y_data[ju].im = y_data[b_i].im - temp_im;
        y_data[b_i].re += temp_re;
        y_data[b_i].im += temp_im;
        b_i += iDelta2;
      }
      istart++;
    }
    k = (int)((unsigned int)k >> 1);
    iDelta = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iDelta;
  }
}

static void d_FFTImplementationCallback_doH(
    const double x_data[], int x_size, emxArray_creal_T *y, int nrowsx,
    int nRows, int nfft, const emxArray_creal_T *wwc,
    const emxArray_real_T *costab, const emxArray_real_T *sintab,
    const emxArray_real_T *costabinv, const emxArray_real_T *sintabinv)
{
  emxArray_creal_T *fv;
  emxArray_creal_T *fy;
  emxArray_creal_T *reconVar1;
  emxArray_creal_T *reconVar2;
  emxArray_int32_T *wrapIndex;
  emxArray_real_T *b_costab;
  emxArray_real_T *b_sintab;
  emxArray_real_T *costab1q;
  emxArray_real_T *hcostabinv;
  emxArray_real_T *hsintab;
  emxArray_real_T *hsintabinv;
  const creal_T *wwc_data;
  creal_T temp;
  creal_T *fv_data;
  creal_T *fy_data;
  creal_T *reconVar1_data;
  creal_T *reconVar2_data;
  creal_T *y_data;
  creal_T *ytmp_data;
  const double *costab_data;
  const double *costabinv_data;
  const double *sintab_data;
  const double *sintabinv_data;
  double b_re_tmp;
  double b_temp_tmp;
  double b_ytmp_im;
  double b_ytmp_re;
  double c_re_tmp;
  double c_ytmp_re_tmp;
  double d;
  double d1;
  double d2;
  double d3;
  double d_re_tmp;
  double d_ytmp_re_tmp;
  double e_ytmp_re_tmp;
  double f_ytmp_re_tmp;
  double g_ytmp_re_tmp;
  double h_ytmp_re_tmp;
  double re_tmp;
  double temp_tmp;
  double twid_im;
  double twid_re;
  double *b_costab_data;
  double *b_sintab_data;
  double *costab1q_data;
  double *hcostabinv_data;
  double *hsintab_data;
  double *hsintabinv_data;
  int b_i;
  int b_k;
  int b_nfft;
  int c_i;
  int c_k;
  int d_i;
  int d_k;
  int e_i;
  int e_k;
  int f_k;
  int hnRows;
  int hszCostab;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int iDelta;
  int iDelta2;
  int iheight;
  int k;
  int n;
  int n2;
  int nRowsD2;
  int nd2;
  int *wrapIndex_data;
  boolean_T tst;
  sintabinv_data = sintabinv->data;
  costabinv_data = costabinv->data;
  sintab_data = sintab->data;
  costab_data = costab->data;
  wwc_data = wwc->data;
  y_data = y->data;
  hnRows = (int)((unsigned int)nRows >> 1);
  emxInit_creal_T(&y);
  nd2 = y->size[0];
  y->size[0] = hnRows;
  emxEnsureCapacity_creal_T(y, nd2);
  ytmp_data = y->data;
  if (hnRows > nrowsx) {
    nd2 = y->size[0];
    y->size[0] = hnRows;
    emxEnsureCapacity_creal_T(y, nd2);
    ytmp_data = y->data;
    for (i = 0; i < hnRows; i++) {
      ytmp_data[i].re = 0.0;
      ytmp_data[i].im = 0.0;
    }
  }
  if (((unsigned int)x_size & 1U) == 0U) {
    tst = true;
    iDelta = x_size;
  } else {
    tst = false;
    iDelta = x_size - 1;
  }
  nd2 = nRows << 1;
  twid_re = 6.2831853071795862 / (double)nd2;
  n = (int)(((unsigned int)nd2 >> 1) >> 1);
  emxInit_real_T(&costab1q, 2);
  nd2 = costab1q->size[0] * costab1q->size[1];
  costab1q->size[0] = 1;
  costab1q->size[1] = n + 1;
  emxEnsureCapacity_real_T(costab1q, nd2);
  costab1q_data = costab1q->data;
  costab1q_data[0] = 1.0;
  nd2 = (int)((unsigned int)n >> 1) - 1;
  if (nd2 + 1 < 1200) {
    for (k = 0; k <= nd2; k++) {
      costab1q_data[k + 1] = cos(twid_re * ((double)k + 1.0));
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (k = 0; k <= nd2; k++) {
      costab1q_data[k + 1] = cos(twid_re * ((double)k + 1.0));
    }
  }
  hszCostab = nd2 + 2;
  if ((n - nd2) - 2 < 1200) {
    for (b_k = hszCostab; b_k < n; b_k++) {
      costab1q_data[b_k] = sin(twid_re * (double)(n - b_k));
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (b_k = hszCostab; b_k < n; b_k++) {
      costab1q_data[b_k] = sin(twid_re * (double)(n - b_k));
    }
  }
  costab1q_data[n] = 0.0;
  n = costab1q->size[1] - 1;
  n2 = (costab1q->size[1] - 1) << 1;
  emxInit_real_T(&b_costab, 2);
  nd2 = b_costab->size[0] * b_costab->size[1];
  b_costab->size[0] = 1;
  b_costab->size[1] = n2 + 1;
  emxEnsureCapacity_real_T(b_costab, nd2);
  b_costab_data = b_costab->data;
  emxInit_real_T(&b_sintab, 2);
  nd2 = b_sintab->size[0] * b_sintab->size[1];
  b_sintab->size[0] = 1;
  b_sintab->size[1] = n2 + 1;
  emxEnsureCapacity_real_T(b_sintab, nd2);
  b_sintab_data = b_sintab->data;
  b_costab_data[0] = 1.0;
  b_sintab_data[0] = 0.0;
  if (costab1q->size[1] - 1 < 1200) {
    for (c_k = 0; c_k < n; c_k++) {
      b_costab_data[c_k + 1] = costab1q_data[c_k + 1];
      b_sintab_data[c_k + 1] = -costab1q_data[(n - c_k) - 1];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (c_k = 0; c_k < n; c_k++) {
      b_costab_data[c_k + 1] = costab1q_data[c_k + 1];
      b_sintab_data[c_k + 1] = -costab1q_data[(n - c_k) - 1];
    }
  }
  hszCostab = costab1q->size[1];
  if ((n2 - costab1q->size[1]) + 1 < 1200) {
    for (d_k = hszCostab; d_k <= n2; d_k++) {
      b_costab_data[d_k] = -costab1q_data[n2 - d_k];
      b_sintab_data[d_k] = -costab1q_data[d_k - n];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (d_k = hszCostab; d_k <= n2; d_k++) {
      b_costab_data[d_k] = -costab1q_data[n2 - d_k];
      b_sintab_data[d_k] = -costab1q_data[d_k - n];
    }
  }
  hszCostab = (int)((unsigned int)costab->size[1] >> 1);
  nd2 = costab1q->size[0] * costab1q->size[1];
  costab1q->size[0] = 1;
  costab1q->size[1] = hszCostab;
  emxEnsureCapacity_real_T(costab1q, nd2);
  costab1q_data = costab1q->data;
  emxInit_real_T(&hsintab, 2);
  nd2 = hsintab->size[0] * hsintab->size[1];
  hsintab->size[0] = 1;
  hsintab->size[1] = hszCostab;
  emxEnsureCapacity_real_T(hsintab, nd2);
  hsintab_data = hsintab->data;
  emxInit_real_T(&hcostabinv, 2);
  nd2 = hcostabinv->size[0] * hcostabinv->size[1];
  hcostabinv->size[0] = 1;
  hcostabinv->size[1] = hszCostab;
  emxEnsureCapacity_real_T(hcostabinv, nd2);
  hcostabinv_data = hcostabinv->data;
  emxInit_real_T(&hsintabinv, 2);
  nd2 = hsintabinv->size[0] * hsintabinv->size[1];
  hsintabinv->size[0] = 1;
  hsintabinv->size[1] = hszCostab;
  emxEnsureCapacity_real_T(hsintabinv, nd2);
  hsintabinv_data = hsintabinv->data;
  if (hszCostab < 1200) {
    for (b_i = 0; b_i < hszCostab; b_i++) {
      c_i = ((b_i + 1) << 1) - 2;
      costab1q_data[b_i] = costab_data[c_i];
      hsintab_data[b_i] = sintab_data[c_i];
      hcostabinv_data[b_i] = costabinv_data[c_i];
      hsintabinv_data[b_i] = sintabinv_data[c_i];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(c_i)

    for (b_i = 0; b_i < hszCostab; b_i++) {
      c_i = ((b_i + 1) << 1) - 2;
      costab1q_data[b_i] = costab_data[c_i];
      hsintab_data[b_i] = sintab_data[c_i];
      hcostabinv_data[b_i] = costabinv_data[c_i];
      hsintabinv_data[b_i] = sintabinv_data[c_i];
    }
  }
  emxInit_creal_T(&reconVar1);
  nd2 = reconVar1->size[0];
  reconVar1->size[0] = hnRows;
  emxEnsureCapacity_creal_T(reconVar1, nd2);
  reconVar1_data = reconVar1->data;
  emxInit_creal_T(&reconVar2);
  nd2 = reconVar2->size[0];
  reconVar2->size[0] = hnRows;
  emxEnsureCapacity_creal_T(reconVar2, nd2);
  reconVar2_data = reconVar2->data;
  emxInit_int32_T(&wrapIndex, 2);
  nd2 = wrapIndex->size[0] * wrapIndex->size[1];
  wrapIndex->size[0] = 1;
  wrapIndex->size[1] = hnRows;
  emxEnsureCapacity_int32_T(wrapIndex, nd2);
  wrapIndex_data = wrapIndex->data;
  for (i = 0; i < hnRows; i++) {
    nd2 = i << 1;
    twid_re = b_sintab_data[nd2];
    twid_im = b_costab_data[nd2];
    reconVar1_data[i].re = twid_re + 1.0;
    reconVar1_data[i].im = -twid_im;
    reconVar2_data[i].re = 1.0 - twid_re;
    reconVar2_data[i].im = twid_im;
    if (i != 0) {
      wrapIndex_data[i] = (hnRows - i) + 1;
    } else {
      wrapIndex_data[0] = 1;
    }
  }
  emxFree_real_T(&b_sintab);
  emxFree_real_T(&b_costab);
  n = (int)((unsigned int)iDelta >> 1);
  nd2 = (unsigned char)n;
  for (i = 0; i < nd2; i++) {
    temp = wwc_data[(hnRows + i) - 1];
    hszCostab = i << 1;
    twid_re = x_data[hszCostab];
    twid_im = x_data[hszCostab + 1];
    ytmp_data[i].re = temp.re * twid_re + temp.im * twid_im;
    ytmp_data[i].im = temp.re * twid_im - temp.im * twid_re;
  }
  if (!tst) {
    temp = wwc_data[(hnRows + n) - 1];
    if ((unsigned char)n - 1 < 0) {
      nd2 = 0;
    } else {
      nd2 = (unsigned char)n << 1;
    }
    twid_re = x_data[nd2];
    ytmp_data[n].re = temp.re * twid_re + temp.im * 0.0;
    ytmp_data[n].im = temp.re * 0.0 - temp.im * twid_re;
    nd2 = n + 2;
    for (i = nd2; i <= hnRows; i++) {
      ytmp_data[i - 1].re = 0.0;
      ytmp_data[i - 1].im = 0.0;
    }
  } else {
    nd2 = n + 1;
    for (i = nd2; i <= hnRows; i++) {
      ytmp_data[i - 1].re = 0.0;
      ytmp_data[i - 1].im = 0.0;
    }
  }
  b_nfft = (int)((unsigned int)nfft >> 1);
  emxInit_creal_T(&fy);
  nd2 = fy->size[0];
  fy->size[0] = b_nfft;
  emxEnsureCapacity_creal_T(fy, nd2);
  fy_data = fy->data;
  if (b_nfft > y->size[0]) {
    nd2 = fy->size[0];
    fy->size[0] = b_nfft;
    emxEnsureCapacity_creal_T(fy, nd2);
    fy_data = fy->data;
    for (i = 0; i < b_nfft; i++) {
      fy_data[i].re = 0.0;
      fy_data[i].im = 0.0;
    }
  }
  n = y->size[0];
  if (n > b_nfft) {
    n = b_nfft;
  }
  iDelta = b_nfft - 2;
  nRowsD2 = (int)((unsigned int)b_nfft >> 1);
  e_k = (int)((unsigned int)nRowsD2 >> 1);
  nd2 = 0;
  n2 = 0;
  for (i = 0; i <= n - 2; i++) {
    fy_data[nd2] = ytmp_data[i];
    nd2 = b_nfft;
    tst = true;
    while (tst) {
      nd2 >>= 1;
      n2 = (int)((unsigned int)n2 ^ (unsigned int)nd2);
      tst = (((unsigned int)n2 & (unsigned int)nd2) == 0U);
    }
    nd2 = n2;
  }
  if (n - 2 < 0) {
    hszCostab = 0;
  } else {
    hszCostab = n - 1;
  }
  fy_data[nd2] = ytmp_data[hszCostab];
  if (b_nfft > 1) {
    for (i = 0; i <= iDelta; i += 2) {
      temp = fy_data[i + 1];
      twid_re = fy_data[i].re;
      twid_im = fy_data[i].im;
      fy_data[i + 1].re = twid_re - fy_data[i + 1].re;
      fy_data[i + 1].im = fy_data[i].im - fy_data[i + 1].im;
      twid_re += temp.re;
      twid_im += temp.im;
      fy_data[i].re = twid_re;
      fy_data[i].im = twid_im;
    }
  }
  iDelta = 2;
  iDelta2 = 4;
  iheight = ((e_k - 1) << 2) + 1;
  while (e_k > 0) {
    for (d_i = 0; d_i < iheight; d_i += iDelta2) {
      nd2 = d_i + iDelta;
      temp = fy_data[nd2];
      fy_data[nd2].re = fy_data[d_i].re - temp.re;
      fy_data[nd2].im = fy_data[d_i].im - temp.im;
      fy_data[d_i].re += temp.re;
      fy_data[d_i].im += temp.im;
    }
    nd2 = 1;
    for (n = e_k; n < nRowsD2; n += e_k) {
      twid_re = costab1q_data[n];
      twid_im = hsintab_data[n];
      d_i = nd2;
      n2 = nd2 + iheight;
      while (d_i < n2) {
        hszCostab = d_i + iDelta;
        temp_tmp = fy_data[hszCostab].im;
        b_temp_tmp = fy_data[hszCostab].re;
        temp.re = twid_re * b_temp_tmp - twid_im * temp_tmp;
        temp.im = twid_re * temp_tmp + twid_im * b_temp_tmp;
        fy_data[hszCostab].re = fy_data[d_i].re - temp.re;
        fy_data[hszCostab].im = fy_data[d_i].im - temp.im;
        fy_data[d_i].re += temp.re;
        fy_data[d_i].im += temp.im;
        d_i += iDelta2;
      }
      nd2++;
    }
    e_k = (int)((unsigned int)e_k >> 1);
    iDelta = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iDelta;
  }
  emxInit_creal_T(&fv);
  c_FFTImplementationCallback_r2b(wwc, b_nfft, costab1q, hsintab, fv);
  fv_data = fv->data;
  emxFree_real_T(&costab1q);
  emxFree_real_T(&hsintab);
  if (fy->size[0] < 1200) {
    for (i1 = 0; i1 < b_nfft; i1++) {
      twid_re = fy_data[i1].re;
      twid_im = fv_data[i1].im;
      temp_tmp = fy_data[i1].im;
      b_temp_tmp = fv_data[i1].re;
      fy_data[i1].re = twid_re * b_temp_tmp - temp_tmp * twid_im;
      fy_data[i1].im = twid_re * twid_im + temp_tmp * b_temp_tmp;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        re_tmp, b_re_tmp, c_re_tmp, d_re_tmp)

    for (i1 = 0; i1 < b_nfft; i1++) {
      re_tmp = fy_data[i1].re;
      b_re_tmp = fv_data[i1].im;
      c_re_tmp = fy_data[i1].im;
      d_re_tmp = fv_data[i1].re;
      fy_data[i1].re = re_tmp * d_re_tmp - c_re_tmp * b_re_tmp;
      fy_data[i1].im = re_tmp * b_re_tmp + c_re_tmp * d_re_tmp;
    }
  }
  c_FFTImplementationCallback_r2b(fy, b_nfft, hcostabinv, hsintabinv, fv);
  fv_data = fv->data;
  emxFree_creal_T(&fy);
  emxFree_real_T(&hsintabinv);
  emxFree_real_T(&hcostabinv);
  if (fv->size[0] > 1) {
    twid_re = 1.0 / (double)fv->size[0];
    nd2 = fv->size[0];
    if (fv->size[0] < 1200) {
      for (i2 = 0; i2 < nd2; i2++) {
        fv_data[i2].re *= twid_re;
        fv_data[i2].im *= twid_re;
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

      for (i2 = 0; i2 < nd2; i2++) {
        fv_data[i2].re *= twid_re;
        fv_data[i2].im *= twid_re;
      }
    }
  }
  hszCostab = wwc->size[0];
  if ((wwc->size[0] - hnRows) + 1 < 1200) {
    for (f_k = hnRows; f_k <= hszCostab; f_k++) {
      twid_re = wwc_data[f_k - 1].re;
      twid_im = fv_data[f_k - 1].im;
      temp_tmp = wwc_data[f_k - 1].im;
      b_temp_tmp = fv_data[f_k - 1].re;
      nd2 = f_k - hnRows;
      ytmp_data[nd2].re = twid_re * b_temp_tmp + temp_tmp * twid_im;
      ytmp_data[nd2].im = twid_re * twid_im - temp_tmp * b_temp_tmp;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        d, d1, d2, d3, i4)

    for (f_k = hnRows; f_k <= hszCostab; f_k++) {
      d = wwc_data[f_k - 1].re;
      d1 = fv_data[f_k - 1].im;
      d2 = wwc_data[f_k - 1].im;
      d3 = fv_data[f_k - 1].re;
      i4 = f_k - hnRows;
      ytmp_data[i4].re = d * d3 + d2 * d1;
      ytmp_data[i4].im = d * d1 - d2 * d3;
    }
  }
  emxFree_creal_T(&fv);
  if (hnRows < 1200) {
    for (e_i = 0; e_i < hnRows; e_i++) {
      double b_ytmp_re_tmp;
      double ytmp_im;
      double ytmp_re;
      double ytmp_re_tmp;
      i3 = wrapIndex_data[e_i];
      twid_re = ytmp_data[e_i].re;
      twid_im = reconVar1_data[e_i].im;
      temp_tmp = ytmp_data[e_i].im;
      b_temp_tmp = reconVar1_data[e_i].re;
      ytmp_re = ytmp_data[i3 - 1].re;
      ytmp_im = -ytmp_data[i3 - 1].im;
      ytmp_re_tmp = reconVar2_data[e_i].im;
      b_ytmp_re_tmp = reconVar2_data[e_i].re;
      y_data[e_i].re =
          0.5 * ((twid_re * b_temp_tmp - temp_tmp * twid_im) +
                 (ytmp_re * b_ytmp_re_tmp - ytmp_im * ytmp_re_tmp));
      y_data[e_i].im =
          0.5 * ((twid_re * twid_im + temp_tmp * b_temp_tmp) +
                 (ytmp_re * ytmp_re_tmp + ytmp_im * b_ytmp_re_tmp));
      nd2 = hnRows + e_i;
      y_data[nd2].re =
          0.5 * ((twid_re * b_ytmp_re_tmp - temp_tmp * ytmp_re_tmp) +
                 (ytmp_re * b_temp_tmp - ytmp_im * twid_im));
      y_data[nd2].im =
          0.5 * ((twid_re * ytmp_re_tmp + temp_tmp * b_ytmp_re_tmp) +
                 (ytmp_re * twid_im + ytmp_im * b_temp_tmp));
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        i3, c_ytmp_re_tmp, d_ytmp_re_tmp, e_ytmp_re_tmp, f_ytmp_re_tmp,        \
            b_ytmp_re, b_ytmp_im, g_ytmp_re_tmp, h_ytmp_re_tmp)

    for (e_i = 0; e_i < hnRows; e_i++) {
      i3 = wrapIndex_data[e_i];
      c_ytmp_re_tmp = ytmp_data[e_i].re;
      d_ytmp_re_tmp = reconVar1_data[e_i].im;
      e_ytmp_re_tmp = ytmp_data[e_i].im;
      f_ytmp_re_tmp = reconVar1_data[e_i].re;
      b_ytmp_re = ytmp_data[i3 - 1].re;
      b_ytmp_im = -ytmp_data[i3 - 1].im;
      g_ytmp_re_tmp = reconVar2_data[e_i].im;
      h_ytmp_re_tmp = reconVar2_data[e_i].re;
      y_data[e_i].re =
          0.5 *
          ((c_ytmp_re_tmp * f_ytmp_re_tmp - e_ytmp_re_tmp * d_ytmp_re_tmp) +
           (b_ytmp_re * h_ytmp_re_tmp - b_ytmp_im * g_ytmp_re_tmp));
      y_data[e_i].im =
          0.5 *
          ((c_ytmp_re_tmp * d_ytmp_re_tmp + e_ytmp_re_tmp * f_ytmp_re_tmp) +
           (b_ytmp_re * g_ytmp_re_tmp + b_ytmp_im * h_ytmp_re_tmp));
      i3 = hnRows + e_i;
      y_data[i3].re =
          0.5 *
          ((c_ytmp_re_tmp * h_ytmp_re_tmp - e_ytmp_re_tmp * g_ytmp_re_tmp) +
           (b_ytmp_re * f_ytmp_re_tmp - b_ytmp_im * d_ytmp_re_tmp));
      y_data[i3].im =
          0.5 *
          ((c_ytmp_re_tmp * g_ytmp_re_tmp + e_ytmp_re_tmp * h_ytmp_re_tmp) +
           (b_ytmp_re * d_ytmp_re_tmp + b_ytmp_im * f_ytmp_re_tmp));
    }
  }
  emxFree_creal_T(&reconVar2);
  emxFree_creal_T(&reconVar1);
  emxFree_int32_T(&wrapIndex);
  emxFree_creal_T(&y);
}

void c_FFTImplementationCallback_doH(const double x_data[], int x_size,
                                     emxArray_creal_T *y, int unsigned_nRows,
                                     const emxArray_real_T *costab,
                                     const emxArray_real_T *sintab)
{
  emxArray_creal_T *reconVar2;
  emxArray_int32_T *bitrevIndex;
  emxArray_int32_T *wrapIndex;
  emxArray_real_T *hcostab;
  emxArray_real_T *hsintab;
  creal_T *reconVar1_data;
  creal_T *reconVar2_data;
  creal_T *y_data;
  const double *costab_data;
  const double *sintab_data;
  double b_im;
  double b_re;
  double im;
  double re;
  double temp2_im;
  double temp2_re;
  double temp_im;
  double temp_re;
  double *hcostab_data;
  double *hsintab_data;
  int b_i;
  int c_i;
  int d_i;
  int hszCostab;
  int i;
  int iDelta;
  int iDelta2;
  int iheight;
  int istart;
  int iy;
  int k;
  int nRows;
  int nRowsD2;
  int u0;
  int *bitrevIndex_data;
  int *wrapIndex_data;
  boolean_T tst;
  sintab_data = sintab->data;
  costab_data = costab->data;
  y_data = y->data;
  nRows = (int)((unsigned int)unsigned_nRows >> 1);
  u0 = y->size[0];
  if (u0 > nRows) {
    u0 = nRows;
  }
  iDelta = u0 - 2;
  iDelta2 = nRows - 2;
  nRowsD2 = (int)((unsigned int)nRows >> 1);
  k = (int)((unsigned int)nRowsD2 >> 1);
  hszCostab = (int)((unsigned int)costab->size[1] >> 1);
  emxInit_real_T(&hcostab, 2);
  istart = hcostab->size[0] * hcostab->size[1];
  hcostab->size[0] = 1;
  hcostab->size[1] = hszCostab;
  emxEnsureCapacity_real_T(hcostab, istart);
  hcostab_data = hcostab->data;
  emxInit_real_T(&hsintab, 2);
  istart = hsintab->size[0] * hsintab->size[1];
  hsintab->size[0] = 1;
  hsintab->size[1] = hszCostab;
  emxEnsureCapacity_real_T(hsintab, istart);
  hsintab_data = hsintab->data;
  if (hszCostab < 1200) {
    for (i = 0; i < hszCostab; i++) {
      b_i = ((i + 1) << 1) - 2;
      hcostab_data[i] = costab_data[b_i];
      hsintab_data[i] = sintab_data[b_i];
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(b_i)

    for (i = 0; i < hszCostab; i++) {
      b_i = ((i + 1) << 1) - 2;
      hcostab_data[i] = costab_data[b_i];
      hsintab_data[i] = sintab_data[b_i];
    }
  }
  emxInit_creal_T(&y);
  istart = y->size[0];
  y->size[0] = nRows;
  emxEnsureCapacity_creal_T(y, istart);
  reconVar1_data = y->data;
  emxInit_creal_T(&reconVar2);
  istart = reconVar2->size[0];
  reconVar2->size[0] = nRows;
  emxEnsureCapacity_creal_T(reconVar2, istart);
  reconVar2_data = reconVar2->data;
  emxInit_int32_T(&wrapIndex, 2);
  istart = wrapIndex->size[0] * wrapIndex->size[1];
  wrapIndex->size[0] = 1;
  wrapIndex->size[1] = nRows;
  emxEnsureCapacity_int32_T(wrapIndex, istart);
  wrapIndex_data = wrapIndex->data;
  hszCostab = 0;
  iy = 1;
  emxInit_int32_T(&bitrevIndex, 1);
  istart = bitrevIndex->size[0];
  bitrevIndex->size[0] = nRows;
  emxEnsureCapacity_int32_T(bitrevIndex, istart);
  bitrevIndex_data = bitrevIndex->data;
  for (c_i = 0; c_i < nRows; c_i++) {
    re = sintab_data[c_i];
    im = costab_data[c_i];
    reconVar1_data[c_i].re = re + 1.0;
    reconVar1_data[c_i].im = -im;
    reconVar2_data[c_i].re = 1.0 - re;
    reconVar2_data[c_i].im = im;
    if (c_i != 0) {
      wrapIndex_data[c_i] = (nRows - c_i) + 1;
    } else {
      wrapIndex_data[0] = 1;
    }
    bitrevIndex_data[c_i] = 0;
  }
  for (c_i = 0; c_i <= iDelta; c_i++) {
    bitrevIndex_data[c_i] = iy;
    istart = nRows;
    tst = true;
    while (tst) {
      istart >>= 1;
      hszCostab = (int)((unsigned int)hszCostab ^ (unsigned int)istart);
      tst = (((unsigned int)hszCostab & (unsigned int)istart) == 0U);
    }
    iy = hszCostab + 1;
  }
  bitrevIndex_data[u0 - 1] = iy;
  if (((unsigned int)x_size & 1U) == 0U) {
    tst = true;
    istart = x_size;
  } else {
    tst = false;
    istart = x_size - 1;
  }
  iy = (int)((unsigned int)istart >> 1);
  istart = (unsigned char)iy;
  for (c_i = 0; c_i < istart; c_i++) {
    hszCostab = c_i << 1;
    y_data[bitrevIndex_data[c_i] - 1].re = x_data[hszCostab];
    y_data[bitrevIndex_data[c_i] - 1].im = x_data[hszCostab + 1];
  }
  if (!tst) {
    if ((unsigned char)iy - 1 < 0) {
      istart = 0;
    } else {
      istart = (unsigned char)iy << 1;
    }
    y_data[bitrevIndex_data[iy] - 1].re = x_data[istart];
    y_data[bitrevIndex_data[iy] - 1].im = 0.0;
  }
  emxFree_int32_T(&bitrevIndex);
  if (nRows > 1) {
    for (c_i = 0; c_i <= iDelta2; c_i += 2) {
      re = y_data[c_i + 1].re;
      im = y_data[c_i + 1].im;
      temp_re = re;
      temp_im = im;
      b_re = y_data[c_i].re;
      b_im = y_data[c_i].im;
      re = b_re - re;
      im = b_im - im;
      y_data[c_i + 1].re = re;
      y_data[c_i + 1].im = im;
      b_re += temp_re;
      b_im += temp_im;
      y_data[c_i].re = b_re;
      y_data[c_i].im = b_im;
    }
  }
  iDelta = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    for (d_i = 0; d_i < iheight; d_i += iDelta2) {
      istart = d_i + iDelta;
      temp_re = y_data[istart].re;
      temp_im = y_data[istart].im;
      y_data[istart].re = y_data[d_i].re - temp_re;
      y_data[istart].im = y_data[d_i].im - temp_im;
      y_data[d_i].re += temp_re;
      y_data[d_i].im += temp_im;
    }
    istart = 1;
    for (hszCostab = k; hszCostab < nRowsD2; hszCostab += k) {
      temp2_re = hcostab_data[hszCostab];
      temp2_im = hsintab_data[hszCostab];
      d_i = istart;
      iy = istart + iheight;
      while (d_i < iy) {
        u0 = d_i + iDelta;
        re = y_data[u0].im;
        im = y_data[u0].re;
        temp_re = temp2_re * im - temp2_im * re;
        temp_im = temp2_re * re + temp2_im * im;
        y_data[u0].re = y_data[d_i].re - temp_re;
        y_data[u0].im = y_data[d_i].im - temp_im;
        y_data[d_i].re += temp_re;
        y_data[d_i].im += temp_im;
        d_i += iDelta2;
      }
      istart++;
    }
    k = (int)((unsigned int)k >> 1);
    iDelta = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iDelta;
  }
  emxFree_real_T(&hsintab);
  emxFree_real_T(&hcostab);
  temp_re = y_data[0].re;
  temp_im = y_data[0].im;
  im = y_data[0].re * reconVar1_data[0].re;
  b_re = y_data[0].re * reconVar1_data[0].im;
  b_im = -y_data[0].im;
  temp2_re = temp_re * reconVar2_data[0].re;
  re = temp_re * reconVar2_data[0].im;
  y_data[0].re = 0.5 * ((im - y_data[0].im * reconVar1_data[0].im) +
                        (temp2_re - b_im * reconVar2_data[0].im));
  y_data[0].im = 0.5 * ((b_re + y_data[0].im * reconVar1_data[0].re) +
                        (re + b_im * reconVar2_data[0].re));
  y_data[nRows].re = 0.5 * ((temp2_re - temp_im * reconVar2_data[0].im) +
                            (im - b_im * reconVar1_data[0].im));
  y_data[nRows].im = 0.5 * ((re + temp_im * reconVar2_data[0].re) +
                            (b_re + b_im * reconVar1_data[0].re));
  for (c_i = 2; c_i <= nRowsD2; c_i++) {
    temp_re = y_data[c_i - 1].re;
    temp_im = y_data[c_i - 1].im;
    istart = wrapIndex_data[c_i - 1];
    temp2_re = y_data[istart - 1].re;
    temp2_im = y_data[istart - 1].im;
    re = reconVar1_data[c_i - 1].im;
    im = reconVar1_data[c_i - 1].re;
    b_re = reconVar2_data[c_i - 1].im;
    b_im = reconVar2_data[c_i - 1].re;
    y_data[c_i - 1].re = 0.5 * ((temp_re * im - temp_im * re) +
                                (temp2_re * b_im - -temp2_im * b_re));
    y_data[c_i - 1].im = 0.5 * ((temp_re * re + temp_im * im) +
                                (temp2_re * b_re + -temp2_im * b_im));
    hszCostab = (nRows + c_i) - 1;
    y_data[hszCostab].re = 0.5 * ((temp_re * b_im - temp_im * b_re) +
                                  (temp2_re * im - -temp2_im * re));
    y_data[hszCostab].im = 0.5 * ((temp_re * b_re + temp_im * b_im) +
                                  (temp2_re * re + -temp2_im * im));
    re = reconVar1_data[istart - 1].im;
    im = reconVar1_data[istart - 1].re;
    b_re = reconVar2_data[istart - 1].im;
    b_im = reconVar2_data[istart - 1].re;
    y_data[istart - 1].re = 0.5 * ((temp2_re * im - temp2_im * re) +
                                   (temp_re * b_im - -temp_im * b_re));
    y_data[istart - 1].im = 0.5 * ((temp2_re * re + temp2_im * im) +
                                   (temp_re * b_re + -temp_im * b_im));
    istart = (istart + nRows) - 1;
    y_data[istart].re = 0.5 * ((temp2_re * b_im - temp2_im * b_re) +
                               (temp_re * im - -temp_im * re));
    y_data[istart].im = 0.5 * ((temp2_re * b_re + temp2_im * b_im) +
                               (temp_re * re + -temp_im * im));
  }
  double b_y_re_tmp;
  double y_re_tmp;
  emxFree_int32_T(&wrapIndex);
  temp_re = y_data[nRowsD2].re;
  temp_im = y_data[nRowsD2].im;
  im = reconVar1_data[nRowsD2].im;
  b_re = reconVar1_data[nRowsD2].re;
  emxFree_creal_T(&y);
  b_im = temp_re * b_re;
  temp2_re = temp_re * im;
  temp2_im = reconVar2_data[nRowsD2].im;
  y_re_tmp = reconVar2_data[nRowsD2].re;
  emxFree_creal_T(&reconVar2);
  b_y_re_tmp = temp_re * y_re_tmp;
  re = temp_re * temp2_im;
  y_data[nRowsD2].re =
      0.5 * ((b_im - temp_im * im) + (b_y_re_tmp - -temp_im * temp2_im));
  y_data[nRowsD2].im =
      0.5 * ((temp2_re + temp_im * b_re) + (re + -temp_im * y_re_tmp));
  istart = nRows + nRowsD2;
  y_data[istart].re =
      0.5 * ((b_y_re_tmp - temp_im * temp2_im) + (b_im - -temp_im * im));
  y_data[istart].im =
      0.5 * ((re + temp_im * y_re_tmp) + (temp2_re + -temp_im * b_re));
}

void c_FFTImplementationCallback_dob(const double x_data[], int x_size,
                                     int n2blue, int nfft,
                                     const emxArray_real_T *costab,
                                     const emxArray_real_T *sintab,
                                     const emxArray_real_T *sintabinv,
                                     emxArray_creal_T *y)
{
  emxArray_creal_T *fv;
  emxArray_creal_T *fy;
  emxArray_creal_T *wwc;
  creal_T *fv_data;
  creal_T *fy_data;
  creal_T *wwc_data;
  creal_T *y_data;
  const double *costab_data;
  const double *sintab_data;
  double b_re_tmp;
  double c_re_tmp;
  double d;
  double d1;
  double d2;
  double d3;
  double d_re_tmp;
  double re_tmp;
  double twid_re;
  int b_i;
  int c_k;
  int i;
  int i1;
  int i2;
  int iDelta;
  int istart;
  int j;
  int k;
  int nInt2;
  int nInt2m1;
  int rt;
  unsigned int u;
  sintab_data = sintab->data;
  costab_data = costab->data;
  u = (unsigned int)nfft & 1U;
  emxInit_creal_T(&wwc);
  if (u == 0U) {
    iDelta = (int)((unsigned int)nfft >> 1);
    nInt2m1 = (iDelta + iDelta) - 1;
    istart = wwc->size[0];
    wwc->size[0] = nInt2m1;
    emxEnsureCapacity_creal_T(wwc, istart);
    wwc_data = wwc->data;
    rt = 0;
    wwc_data[iDelta - 1].re = 1.0;
    wwc_data[iDelta - 1].im = 0.0;
    nInt2 = iDelta << 1;
    for (k = 0; k <= iDelta - 2; k++) {
      istart = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= istart) {
        rt += istart - nInt2;
      } else {
        rt += istart;
      }
      twid_re = -3.1415926535897931 * (double)rt / (double)iDelta;
      j = (iDelta - k) - 2;
      wwc_data[j].re = cos(twid_re);
      wwc_data[j].im = -sin(twid_re);
    }
    istart = nInt2m1 - 1;
    for (k = istart; k >= iDelta; k--) {
      wwc_data[k] = wwc_data[(nInt2m1 - k) - 1];
    }
  } else {
    nInt2 = (nfft + nfft) - 1;
    istart = wwc->size[0];
    wwc->size[0] = nInt2;
    emxEnsureCapacity_creal_T(wwc, istart);
    wwc_data = wwc->data;
    j = 0;
    wwc_data[nfft - 1].re = 1.0;
    wwc_data[nfft - 1].im = 0.0;
    rt = nfft << 1;
    for (k = 0; k <= nfft - 2; k++) {
      istart = ((k + 1) << 1) - 1;
      if (rt - j <= istart) {
        j += istart - rt;
      } else {
        j += istart;
      }
      twid_re = -3.1415926535897931 * (double)j / (double)nfft;
      istart = (nfft - k) - 2;
      wwc_data[istart].re = cos(twid_re);
      wwc_data[istart].im = -sin(twid_re);
    }
    istart = nInt2 - 1;
    for (k = istart; k >= nfft; k--) {
      wwc_data[k] = wwc_data[(nInt2 - k) - 1];
    }
  }
  istart = y->size[0];
  y->size[0] = nfft;
  emxEnsureCapacity_creal_T(y, istart);
  y_data = y->data;
  if (nfft > x_size) {
    istart = y->size[0];
    y->size[0] = nfft;
    emxEnsureCapacity_creal_T(y, istart);
    y_data = y->data;
    for (k = 0; k < nfft; k++) {
      y_data[k].re = 0.0;
      y_data[k].im = 0.0;
    }
  }
  emxInit_creal_T(&fy);
  emxInit_creal_T(&fv);
  if ((n2blue != 1) && (u == 0U)) {
    d_FFTImplementationCallback_doH(x_data, x_size, y, x_size, nfft, n2blue,
                                    wwc, costab, sintab, costab, sintabinv);
  } else {
    creal_T nt;
    double b_nt_tmp;
    double nt_tmp;
    double twid_im;
    int b_k;
    int iheight;
    int nRowsD2;
    for (k = 0; k < x_size; k++) {
      nt = wwc_data[(nfft + k) - 1];
      twid_re = x_data[k];
      y_data[k].re = nt.re * twid_re;
      y_data[k].im = nt.im * -twid_re;
    }
    istart = x_size + 1;
    for (k = istart; k <= nfft; k++) {
      y_data[k - 1].re = 0.0;
      y_data[k - 1].im = 0.0;
    }
    istart = fy->size[0];
    fy->size[0] = n2blue;
    emxEnsureCapacity_creal_T(fy, istart);
    fy_data = fy->data;
    if (n2blue > y->size[0]) {
      istart = fy->size[0];
      fy->size[0] = n2blue;
      emxEnsureCapacity_creal_T(fy, istart);
      fy_data = fy->data;
      for (k = 0; k < n2blue; k++) {
        fy_data[k].re = 0.0;
        fy_data[k].im = 0.0;
      }
    }
    j = y->size[0];
    if (j > n2blue) {
      j = n2blue;
    }
    iDelta = n2blue - 2;
    nRowsD2 = (int)((unsigned int)n2blue >> 1);
    b_k = (int)((unsigned int)nRowsD2 >> 1);
    rt = 0;
    nInt2 = 0;
    for (k = 0; k <= j - 2; k++) {
      boolean_T tst;
      fy_data[rt] = y_data[k];
      istart = n2blue;
      tst = true;
      while (tst) {
        istart >>= 1;
        nInt2 ^= istart;
        tst = ((nInt2 & istart) == 0);
      }
      rt = nInt2;
    }
    if (j - 2 < 0) {
      istart = 0;
    } else {
      istart = j - 1;
    }
    fy_data[rt] = y_data[istart];
    if (n2blue > 1) {
      for (k = 0; k <= iDelta; k += 2) {
        nt = fy_data[k + 1];
        twid_re = fy_data[k].re;
        twid_im = fy_data[k].im;
        fy_data[k + 1].re = twid_re - fy_data[k + 1].re;
        fy_data[k + 1].im = fy_data[k].im - fy_data[k + 1].im;
        twid_re += nt.re;
        twid_im += nt.im;
        fy_data[k].re = twid_re;
        fy_data[k].im = twid_im;
      }
    }
    iDelta = 2;
    nInt2m1 = 4;
    iheight = ((b_k - 1) << 2) + 1;
    while (b_k > 0) {
      for (i = 0; i < iheight; i += nInt2m1) {
        istart = i + iDelta;
        nt = fy_data[istart];
        fy_data[istart].re = fy_data[i].re - nt.re;
        fy_data[istart].im = fy_data[i].im - nt.im;
        fy_data[i].re += nt.re;
        fy_data[i].im += nt.im;
      }
      istart = 1;
      for (j = b_k; j < nRowsD2; j += b_k) {
        twid_re = costab_data[j];
        twid_im = sintab_data[j];
        i = istart;
        rt = istart + iheight;
        while (i < rt) {
          nInt2 = i + iDelta;
          nt_tmp = fy_data[nInt2].im;
          b_nt_tmp = fy_data[nInt2].re;
          nt.re = twid_re * b_nt_tmp - twid_im * nt_tmp;
          nt.im = twid_re * nt_tmp + twid_im * b_nt_tmp;
          fy_data[nInt2].re = fy_data[i].re - nt.re;
          fy_data[nInt2].im = fy_data[i].im - nt.im;
          fy_data[i].re += nt.re;
          fy_data[i].im += nt.im;
          i += nInt2m1;
        }
        istart++;
      }
      b_k = (int)((unsigned int)b_k >> 1);
      iDelta = nInt2m1;
      nInt2m1 += nInt2m1;
      iheight -= iDelta;
    }
    c_FFTImplementationCallback_r2b(wwc, n2blue, costab, sintab, fv);
    fv_data = fv->data;
    istart = fy->size[0];
    if (fy->size[0] < 1200) {
      for (b_i = 0; b_i < istart; b_i++) {
        twid_re = fy_data[b_i].re;
        twid_im = fv_data[b_i].im;
        nt_tmp = fy_data[b_i].im;
        b_nt_tmp = fv_data[b_i].re;
        fy_data[b_i].re = twid_re * b_nt_tmp - nt_tmp * twid_im;
        fy_data[b_i].im = twid_re * twid_im + nt_tmp * b_nt_tmp;
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        re_tmp, b_re_tmp, c_re_tmp, d_re_tmp)

      for (b_i = 0; b_i < istart; b_i++) {
        re_tmp = fy_data[b_i].re;
        b_re_tmp = fv_data[b_i].im;
        c_re_tmp = fy_data[b_i].im;
        d_re_tmp = fv_data[b_i].re;
        fy_data[b_i].re = re_tmp * d_re_tmp - c_re_tmp * b_re_tmp;
        fy_data[b_i].im = re_tmp * b_re_tmp + c_re_tmp * d_re_tmp;
      }
    }
    c_FFTImplementationCallback_r2b(fy, n2blue, costab, sintabinv, fv);
    fv_data = fv->data;
    if (fv->size[0] > 1) {
      twid_re = 1.0 / (double)fv->size[0];
      istart = fv->size[0];
      if (fv->size[0] < 1200) {
        for (i1 = 0; i1 < istart; i1++) {
          fv_data[i1].re *= twid_re;
          fv_data[i1].im *= twid_re;
        }
      } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

        for (i1 = 0; i1 < istart; i1++) {
          fv_data[i1].re *= twid_re;
          fv_data[i1].im *= twid_re;
        }
      }
    }
    j = wwc->size[0];
    if ((wwc->size[0] - nfft) + 1 < 1200) {
      for (c_k = nfft; c_k <= j; c_k++) {
        twid_re = wwc_data[c_k - 1].re;
        twid_im = fv_data[c_k - 1].im;
        nt_tmp = wwc_data[c_k - 1].im;
        b_nt_tmp = fv_data[c_k - 1].re;
        istart = c_k - nfft;
        y_data[istart].re = twid_re * b_nt_tmp + nt_tmp * twid_im;
        y_data[istart].im = twid_re * twid_im - nt_tmp * b_nt_tmp;
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        d, d1, d2, d3, i2)

      for (c_k = nfft; c_k <= j; c_k++) {
        d = wwc_data[c_k - 1].re;
        d1 = fv_data[c_k - 1].im;
        d2 = wwc_data[c_k - 1].im;
        d3 = fv_data[c_k - 1].re;
        i2 = c_k - nfft;
        y_data[i2].re = d * d3 + d2 * d1;
        y_data[i2].im = d * d1 - d2 * d3;
      }
    }
  }
  emxFree_creal_T(&fv);
  emxFree_creal_T(&fy);
  emxFree_creal_T(&wwc);
}

/* End of code generation (FFTImplementationCallback.c) */
