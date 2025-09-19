/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * pca.c
 *
 * Code generation for function 'pca'
 *
 */

/* Include files */
#include "pca.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "xzsvdc.h"
#include "omp.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Definitions */
int pca(const emxArray_real_T *x, double varargout_1_data[],
        int varargout_1_size[2], emxArray_real_T *varargout_2,
        double varargout_3_data[])
{
  __m128d r;
  __m128d r1;
  emxArray_boolean_T *naninfo_isNaN;
  emxArray_int32_T *naninfo_nNaNsInRow;
  emxArray_real_T *b_x;
  emxArray_real_T *score;
  emxArray_real_T *xNoNaNs;
  emxArray_real_T *y;
  double b_coeff_data[9];
  double coeff_data[9];
  double latent_data[3];
  double mu[3];
  const double *x_data;
  double d;
  double maxval;
  double wcol;
  double xcol;
  double *score_data;
  double *xNoNaNs_data;
  double *y_data;
  int coeff_size[2];
  int DOF;
  int b_i;
  int b_j;
  int b_n;
  int c_i;
  int c_j;
  int c_n;
  int d_i;
  int d_j;
  int e_i;
  int i;
  int irow;
  int j;
  int m;
  int n;
  int nrows;
  int varargout_3_size;
  int *naninfo_nNaNsInRow_data;
  boolean_T noNaNs;
  boolean_T *naninfo_isNaN_data;
  x_data = x->data;
  emxInit_real_T(&b_x, 2);
  m = x->size[0];
  irow = b_x->size[0] * b_x->size[1];
  b_x->size[0] = x->size[0];
  b_x->size[1] = 3;
  emxEnsureCapacity_real_T(b_x, irow);
  y_data = b_x->data;
  n = x->size[0] * 3;
  for (i = 0; i < n; i++) {
    y_data[i] = x_data[i];
  }
  b_n = x->size[0];
  c_n = x->size[0];
  nrows = 0;
  varargout_3_size = 0;
  emxInit_int32_T(&naninfo_nNaNsInRow, 1);
  irow = naninfo_nNaNsInRow->size[0];
  naninfo_nNaNsInRow->size[0] = x->size[0];
  emxEnsureCapacity_int32_T(naninfo_nNaNsInRow, irow);
  naninfo_nNaNsInRow_data = naninfo_nNaNsInRow->data;
  for (i = 0; i < m; i++) {
    naninfo_nNaNsInRow_data[i] = 0;
  }
  emxInit_boolean_T(&naninfo_isNaN, 2);
  irow = naninfo_isNaN->size[0] * naninfo_isNaN->size[1];
  naninfo_isNaN->size[0] = x->size[0];
  naninfo_isNaN->size[1] = 3;
  emxEnsureCapacity_boolean_T(naninfo_isNaN, irow);
  naninfo_isNaN_data = naninfo_isNaN->data;
  irow = (n < 1200);
  if (irow) {
    for (b_i = 0; b_i < n; b_i++) {
      naninfo_isNaN_data[b_i] = rtIsNaN(x_data[b_i]);
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads())

    for (b_i = 0; b_i < n; b_i++) {
      naninfo_isNaN_data[b_i] = rtIsNaN(x_data[b_i]);
    }
  }
  for (j = 0; j < 3; j++) {
    for (i = 0; i < c_n; i++) {
      if (naninfo_isNaN_data[i + naninfo_isNaN->size[0] * j]) {
        naninfo_nNaNsInRow_data[i]++;
        nrows++;
      }
    }
  }
  emxFree_boolean_T(&naninfo_isNaN);
  for (i = 0; i < c_n; i++) {
    if (naninfo_nNaNsInRow_data[i] > 0) {
      varargout_3_size++;
    }
  }
  noNaNs = (nrows <= 0);
  DOF = x->size[0] - varargout_3_size;
  if (DOF >= 1) {
    DOF--;
  }
  if (!noNaNs) {
    if (irow) {
      for (c_j = 0; c_j < 3; c_j++) {
        wcol = 0.0;
        xcol = 0.0;
        for (d_i = 0; d_i < m; d_i++) {
          d = x_data[d_i + x->size[0] * c_j];
          if (!rtIsNaN(d)) {
            wcol++;
            xcol += d;
          }
        }
        mu[c_j] = xcol / wcol;
      }
    } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        d, xcol, wcol, d_i)

      for (c_j = 0; c_j < 3; c_j++) {
        wcol = 0.0;
        xcol = 0.0;
        for (d_i = 0; d_i < m; d_i++) {
          d = x_data[d_i + x->size[0] * c_j];
          if (!rtIsNaN(d)) {
            wcol++;
            xcol += d;
          }
        }
        mu[c_j] = xcol / wcol;
      }
    }
  } else if (irow) {
    for (b_j = 0; b_j < 3; b_j++) {
      wcol = 0.0;
      xcol = 0.0;
      for (c_i = 0; c_i < m; c_i++) {
        wcol++;
        xcol += x_data[c_i + x->size[0] * b_j];
      }
      mu[b_j] = xcol / wcol;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        xcol, wcol, c_i)

    for (b_j = 0; b_j < 3; b_j++) {
      wcol = 0.0;
      xcol = 0.0;
      for (c_i = 0; c_i < m; c_i++) {
        wcol++;
        xcol += x_data[c_i + x->size[0] * b_j];
      }
      mu[b_j] = xcol / wcol;
    }
  }
  c_n = (b_n / 2) << 1;
  nrows = c_n - 2;
  irow = (c_n - 1) / 2;
  n = b_n - c_n;
  if (irow >= n) {
    n = irow;
  }
  if (3 * n < 1200) {
    for (d_j = 0; d_j < 3; d_j++) {
      for (e_i = 0; e_i <= nrows; e_i += 2) {
        r = _mm_loadu_pd(&y_data[e_i + b_x->size[0] * d_j]);
        _mm_storeu_pd(&y_data[e_i + b_x->size[0] * d_j],
                      _mm_sub_pd(r, _mm_set1_pd(mu[d_j])));
      }
      for (e_i = c_n; e_i < b_n; e_i++) {
        y_data[e_i + b_x->size[0] * d_j] -= mu[d_j];
      }
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(r, e_i)

    for (d_j = 0; d_j < 3; d_j++) {
      for (e_i = 0; e_i <= nrows; e_i += 2) {
        r = _mm_loadu_pd(&y_data[e_i + b_x->size[0] * d_j]);
        r = _mm_sub_pd(r, _mm_set1_pd(mu[d_j]));
        _mm_storeu_pd(&y_data[e_i + b_x->size[0] * d_j], r);
      }
      for (e_i = c_n; e_i < b_n; e_i++) {
        y_data[e_i + b_x->size[0] * d_j] -= mu[d_j];
      }
    }
  }
  emxInit_real_T(&xNoNaNs, 2);
  emxInit_real_T(&y, 2);
  emxInit_real_T(&score, 2);
  if (noNaNs) {
    nrows = b_x->size[0];
    irow = xNoNaNs->size[0] * xNoNaNs->size[1];
    xNoNaNs->size[0] = b_x->size[0];
    xNoNaNs->size[1] = 3;
    emxEnsureCapacity_real_T(xNoNaNs, irow);
    xNoNaNs_data = xNoNaNs->data;
    irow = b_x->size[0] * b_x->size[1] - 1;
    for (i = 0; i <= irow; i++) {
      xNoNaNs_data[i] = y_data[i];
    }
    varargout_3_size = xzsvdc(xNoNaNs, y, latent_data, coeff_data, coeff_size);
    y_data = y->data;
    c_n = y->size[1];
    for (j = 0; j < c_n; j++) {
      irow = (nrows / 2) << 1;
      n = irow - 2;
      for (i = 0; i <= n; i += 2) {
        r1 = _mm_loadu_pd(&y_data[i + y->size[0] * j]);
        _mm_storeu_pd(&y_data[i + y->size[0] * j],
                      _mm_mul_pd(r1, _mm_set1_pd(latent_data[j])));
      }
      for (i = irow; i < nrows; i++) {
        y_data[i + y->size[0] * j] *= latent_data[j];
      }
    }
    irow = (y->size[1] / 2) << 1;
    n = irow - 2;
    for (j = 0; j <= n; j += 2) {
      r1 = _mm_loadu_pd(&latent_data[0]);
      _mm_storeu_pd(&latent_data[0],
                    _mm_div_pd(_mm_mul_pd(r1, r1), _mm_set1_pd(DOF)));
    }
    for (j = irow; j < c_n; j++) {
      maxval = latent_data[j];
      maxval = maxval * maxval / (double)DOF;
      latent_data[j] = maxval;
    }
    if (DOF < 3) {
      b_n = y->size[1];
      if (DOF <= b_n) {
        b_n = DOF;
      }
      irow = score->size[0] * score->size[1];
      score->size[0] = b_x->size[0];
      score->size[1] = b_n;
      emxEnsureCapacity_real_T(score, irow);
      score_data = score->data;
      for (i = 0; i < b_n; i++) {
        for (j = 0; j < nrows; j++) {
          score_data[j + score->size[0] * i] = y_data[j + y->size[0] * i];
        }
      }
      varargout_3_size = b_n;
      for (j = 0; j < b_n; j++) {
        varargout_3_data[j] = latent_data[j];
        b_coeff_data[3 * j] = coeff_data[3 * j];
        irow = 3 * j + 1;
        b_coeff_data[irow] = coeff_data[irow];
        irow = 3 * j + 2;
        b_coeff_data[irow] = coeff_data[irow];
      }
    } else {
      irow = score->size[0] * score->size[1];
      score->size[0] = y->size[0];
      score->size[1] = y->size[1];
      emxEnsureCapacity_real_T(score, irow);
      score_data = score->data;
      irow = y->size[0] * y->size[1];
      for (i = 0; i < irow; i++) {
        score_data[i] = y_data[i];
      }
      if (varargout_3_size - 1 >= 0) {
        memcpy(&varargout_3_data[0], &latent_data[0],
               (unsigned int)varargout_3_size * sizeof(double));
      }
      b_n = coeff_size[1];
      irow = 3 * coeff_size[1];
      if (irow - 1 >= 0) {
        memcpy(&b_coeff_data[0], &coeff_data[0],
               (unsigned int)irow * sizeof(double));
      }
    }
  } else {
    n = b_x->size[0];
    c_n = b_x->size[0] - varargout_3_size;
    irow = xNoNaNs->size[0] * xNoNaNs->size[1];
    xNoNaNs->size[0] = c_n;
    xNoNaNs->size[1] = 3;
    emxEnsureCapacity_real_T(xNoNaNs, irow);
    xNoNaNs_data = xNoNaNs->data;
    irow = -1;
    for (i = 0; i < n; i++) {
      if (naninfo_nNaNsInRow_data[i] == 0) {
        irow++;
        xNoNaNs_data[irow] = y_data[i];
        xNoNaNs_data[irow + xNoNaNs->size[0]] = y_data[i + b_x->size[0]];
        xNoNaNs_data[irow + xNoNaNs->size[0] * 2] =
            y_data[i + b_x->size[0] * 2];
      }
    }
    varargout_3_size =
        xzsvdc(xNoNaNs, score, latent_data, coeff_data, coeff_size);
    score_data = score->data;
    nrows = score->size[1];
    for (j = 0; j < nrows; j++) {
      irow = (c_n / 2) << 1;
      n = irow - 2;
      for (i = 0; i <= n; i += 2) {
        r1 = _mm_loadu_pd(&score_data[i + score->size[0] * j]);
        _mm_storeu_pd(&score_data[i + score->size[0] * j],
                      _mm_mul_pd(r1, _mm_set1_pd(latent_data[j])));
      }
      for (i = irow; i < c_n; i++) {
        score_data[i + score->size[0] * j] *= latent_data[j];
      }
    }
    irow = (score->size[1] / 2) << 1;
    n = irow - 2;
    for (j = 0; j <= n; j += 2) {
      r1 = _mm_loadu_pd(&latent_data[0]);
      _mm_storeu_pd(&latent_data[0],
                    _mm_div_pd(_mm_mul_pd(r1, r1), _mm_set1_pd(DOF)));
    }
    for (j = irow; j < nrows; j++) {
      maxval = latent_data[j];
      maxval = maxval * maxval / (double)DOF;
      latent_data[j] = maxval;
    }
    if (DOF < 3) {
      b_n = score->size[1];
      if (DOF <= b_n) {
        b_n = DOF;
      }
      irow = y->size[0] * y->size[1];
      y->size[0] = c_n;
      y->size[1] = b_n;
      emxEnsureCapacity_real_T(y, irow);
      y_data = y->data;
      for (i = 0; i < b_n; i++) {
        for (j = 0; j < c_n; j++) {
          y_data[j + y->size[0] * i] = score_data[j + score->size[0] * i];
        }
      }
      varargout_3_size = b_n;
      for (j = 0; j < b_n; j++) {
        varargout_3_data[j] = latent_data[j];
        b_coeff_data[3 * j] = coeff_data[3 * j];
        irow = 3 * j + 1;
        b_coeff_data[irow] = coeff_data[irow];
        irow = 3 * j + 2;
        b_coeff_data[irow] = coeff_data[irow];
      }
    } else {
      irow = y->size[0] * y->size[1];
      y->size[0] = score->size[0];
      y->size[1] = score->size[1];
      emxEnsureCapacity_real_T(y, irow);
      y_data = y->data;
      irow = score->size[0] * score->size[1];
      for (i = 0; i < irow; i++) {
        y_data[i] = score_data[i];
      }
      if (varargout_3_size - 1 >= 0) {
        memcpy(&varargout_3_data[0], &latent_data[0],
               (unsigned int)varargout_3_size * sizeof(double));
      }
      b_n = coeff_size[1];
      irow = 3 * coeff_size[1];
      if (irow - 1 >= 0) {
        memcpy(&b_coeff_data[0], &coeff_data[0],
               (unsigned int)irow * sizeof(double));
      }
    }
    n = y->size[1];
    irow = score->size[0] * score->size[1];
    score->size[0] = m;
    score->size[1] = y->size[1];
    emxEnsureCapacity_real_T(score, irow);
    score_data = score->data;
    irow = -1;
    for (i = 0; i < m; i++) {
      if (naninfo_nNaNsInRow_data[i] > 0) {
        for (j = 0; j < n; j++) {
          score_data[i + score->size[0] * j] = rtNaN;
        }
      } else {
        irow++;
        for (j = 0; j < n; j++) {
          score_data[i + score->size[0] * j] = y_data[irow + y->size[0] * j];
        }
      }
    }
  }
  emxFree_real_T(&y);
  emxFree_real_T(&xNoNaNs);
  emxFree_int32_T(&naninfo_nNaNsInRow);
  emxFree_real_T(&b_x);
  nrows = score->size[0];
  if (DOF > 3) {
    varargout_1_size[0] = 3;
    varargout_1_size[1] = 3;
    for (i = 0; i < 3; i++) {
      varargout_1_data[3 * i] = b_coeff_data[3 * i];
      irow = 3 * i + 1;
      varargout_1_data[irow] = b_coeff_data[irow];
      irow = 3 * i + 2;
      varargout_1_data[irow] = b_coeff_data[irow];
    }
    irow = varargout_2->size[0] * varargout_2->size[1];
    varargout_2->size[0] = score->size[0];
    varargout_2->size[1] = 3;
    emxEnsureCapacity_real_T(varargout_2, irow);
    xNoNaNs_data = varargout_2->data;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < nrows; j++) {
        xNoNaNs_data[j + varargout_2->size[0] * i] =
            score_data[j + score->size[0] * i];
      }
    }
  } else {
    varargout_1_size[0] = 3;
    varargout_1_size[1] = b_n;
    irow = 3 * b_n;
    if (irow - 1 >= 0) {
      memcpy(&varargout_1_data[0], &b_coeff_data[0],
             (unsigned int)irow * sizeof(double));
    }
    irow = varargout_2->size[0] * varargout_2->size[1];
    varargout_2->size[0] = score->size[0];
    varargout_2->size[1] = score->size[1];
    emxEnsureCapacity_real_T(varargout_2, irow);
    xNoNaNs_data = varargout_2->data;
    irow = score->size[0] * score->size[1];
    for (i = 0; i < irow; i++) {
      xNoNaNs_data[i] = score_data[i];
    }
  }
  emxFree_real_T(&score);
  c_n = varargout_1_size[1];
  for (i = 0; i < c_n; i++) {
    double absc;
    double d1;
    double sgn;
    maxval = 0.0;
    sgn = 1.0;
    d1 = varargout_1_data[3 * i];
    absc = fabs(d1);
    if (absc > 0.0) {
      maxval = absc;
      sgn = d1;
      if (!rtIsNaN(d1)) {
        if (d1 < 0.0) {
          sgn = -1.0;
        } else {
          sgn = (d1 > 0.0);
        }
      }
    }
    d1 = varargout_1_data[3 * i + 1];
    absc = fabs(d1);
    if (absc > maxval) {
      maxval = absc;
      sgn = d1;
      if (!rtIsNaN(d1)) {
        if (d1 < 0.0) {
          sgn = -1.0;
        } else {
          sgn = (d1 > 0.0);
        }
      }
    }
    irow = 3 * i + 2;
    d1 = varargout_1_data[irow];
    if (fabs(d1) > maxval) {
      sgn = d1;
      if (!rtIsNaN(d1)) {
        if (d1 < 0.0) {
          sgn = -1.0;
        } else {
          sgn = (d1 > 0.0);
        }
      }
    }
    if (sgn < 0.0) {
      __m128d r2;
      r1 = _mm_loadu_pd(&varargout_1_data[3 * i]);
      r2 = _mm_set1_pd(-1.0);
      _mm_storeu_pd(&varargout_1_data[3 * i], _mm_mul_pd(r1, r2));
      varargout_1_data[irow] = -varargout_1_data[irow];
      irow = (nrows / 2) << 1;
      n = irow - 2;
      for (j = 0; j <= n; j += 2) {
        r1 = _mm_loadu_pd(&xNoNaNs_data[j + varargout_2->size[0] * i]);
        _mm_storeu_pd(&xNoNaNs_data[j + varargout_2->size[0] * i],
                      _mm_mul_pd(r1, r2));
      }
      for (j = irow; j < nrows; j++) {
        xNoNaNs_data[j + varargout_2->size[0] * i] =
            -xNoNaNs_data[j + varargout_2->size[0] * i];
      }
    }
  }
  return varargout_3_size;
}

/* End of code generation (pca.c) */
