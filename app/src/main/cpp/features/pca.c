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
  double wcol;
  double xcol;
  double *score_data;
  double *xNoNaNs_data;
  double *y_data;
  int coeff_size[2];
  int DOF;
  int i;
  int irow;
  int j;
  int loop_ub;
  int n;
  int nrows;
  int nsv;
  int varargout_3_size;
  int *naninfo_nNaNsInRow_data;
  boolean_T noNaNs;
  boolean_T *naninfo_isNaN_data;
  x_data = x->data;
  emxInit_real_T(&b_x, 2);
  loop_ub = x->size[0];
  irow = b_x->size[0] * b_x->size[1];
  b_x->size[0] = x->size[0];
  b_x->size[1] = 3;
  emxEnsureCapacity_real_T(b_x, irow);
  y_data = b_x->data;
  n = x->size[0] * 3;
  for (j = 0; j < n; j++) {
    y_data[j] = x_data[j];
  }
  nsv = 0;
  nrows = 0;
  emxInit_int32_T(&naninfo_nNaNsInRow, 1);
  irow = naninfo_nNaNsInRow->size[0];
  naninfo_nNaNsInRow->size[0] = x->size[0];
  emxEnsureCapacity_int32_T(naninfo_nNaNsInRow, irow);
  naninfo_nNaNsInRow_data = naninfo_nNaNsInRow->data;
  for (j = 0; j < loop_ub; j++) {
    naninfo_nNaNsInRow_data[j] = 0;
  }
  emxInit_boolean_T(&naninfo_isNaN, 2);
  irow = naninfo_isNaN->size[0] * naninfo_isNaN->size[1];
  naninfo_isNaN->size[0] = x->size[0];
  naninfo_isNaN->size[1] = 3;
  emxEnsureCapacity_boolean_T(naninfo_isNaN, irow);
  naninfo_isNaN_data = naninfo_isNaN->data;
  for (j = 0; j < n; j++) {
    naninfo_isNaN_data[j] = rtIsNaN(x_data[j]);
  }
  for (j = 0; j < 3; j++) {
    for (i = 0; i < loop_ub; i++) {
      if (naninfo_isNaN_data[i + naninfo_isNaN->size[0] * j]) {
        naninfo_nNaNsInRow_data[i]++;
        nsv++;
      }
    }
  }
  emxFree_boolean_T(&naninfo_isNaN);
  for (j = 0; j < loop_ub; j++) {
    if (naninfo_nNaNsInRow_data[j] > 0) {
      nrows++;
    }
  }
  noNaNs = (nsv <= 0);
  DOF = x->size[0] - nrows;
  if (DOF >= 1) {
    DOF--;
  }
  if (!noNaNs) {
    for (j = 0; j < 3; j++) {
      wcol = 0.0;
      xcol = 0.0;
      for (i = 0; i < loop_ub; i++) {
        d = x_data[i + x->size[0] * j];
        if (!rtIsNaN(d)) {
          wcol++;
          xcol += d;
        }
      }
      mu[j] = xcol / wcol;
    }
  } else {
    for (j = 0; j < 3; j++) {
      wcol = 0.0;
      xcol = 0.0;
      for (i = 0; i < loop_ub; i++) {
        wcol++;
        xcol += x_data[i + x->size[0] * j];
      }
      mu[j] = xcol / wcol;
    }
  }
  irow = (x->size[0] / 2) << 1;
  n = irow - 2;
  for (j = 0; j < 3; j++) {
    for (i = 0; i <= n; i += 2) {
      r = _mm_loadu_pd(&y_data[i + b_x->size[0] * j]);
      _mm_storeu_pd(&y_data[i + b_x->size[0] * j],
                    _mm_sub_pd(r, _mm_set1_pd(mu[j])));
    }
    for (i = irow; i < loop_ub; i++) {
      y_data[i + b_x->size[0] * j] -= mu[j];
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
    for (j = 0; j <= irow; j++) {
      xNoNaNs_data[j] = y_data[j];
    }
    varargout_3_size = xzsvdc(xNoNaNs, y, latent_data, coeff_data, coeff_size);
    y_data = y->data;
    nsv = y->size[1];
    for (j = 0; j < nsv; j++) {
      irow = (nrows / 2) << 1;
      n = irow - 2;
      for (i = 0; i <= n; i += 2) {
        r = _mm_loadu_pd(&y_data[i + y->size[0] * j]);
        _mm_storeu_pd(&y_data[i + y->size[0] * j],
                      _mm_mul_pd(r, _mm_set1_pd(latent_data[j])));
      }
      for (i = irow; i < nrows; i++) {
        y_data[i + y->size[0] * j] *= latent_data[j];
      }
    }
    irow = (y->size[1] / 2) << 1;
    n = irow - 2;
    for (j = 0; j <= n; j += 2) {
      r = _mm_loadu_pd(&latent_data[0]);
      _mm_storeu_pd(&latent_data[0],
                    _mm_div_pd(_mm_mul_pd(r, r), _mm_set1_pd(DOF)));
    }
    for (j = irow; j < nsv; j++) {
      wcol = latent_data[j];
      wcol = wcol * wcol / (double)DOF;
      latent_data[j] = wcol;
    }
    if (DOF < 3) {
      nsv = y->size[1];
      if (DOF <= nsv) {
        nsv = DOF;
      }
      irow = score->size[0] * score->size[1];
      score->size[0] = b_x->size[0];
      score->size[1] = nsv;
      emxEnsureCapacity_real_T(score, irow);
      score_data = score->data;
      for (j = 0; j < nsv; j++) {
        for (i = 0; i < nrows; i++) {
          score_data[i + score->size[0] * j] = y_data[i + y->size[0] * j];
        }
      }
      varargout_3_size = nsv;
      for (j = 0; j < nsv; j++) {
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
      for (j = 0; j < irow; j++) {
        score_data[j] = y_data[j];
      }
      if (varargout_3_size - 1 >= 0) {
        memcpy(&varargout_3_data[0], &latent_data[0],
               (unsigned int)varargout_3_size * sizeof(double));
      }
      nsv = coeff_size[1];
      irow = 3 * coeff_size[1];
      if (irow - 1 >= 0) {
        memcpy(&b_coeff_data[0], &coeff_data[0],
               (unsigned int)irow * sizeof(double));
      }
    }
  } else {
    n = b_x->size[0];
    nrows = b_x->size[0] - nrows;
    irow = xNoNaNs->size[0] * xNoNaNs->size[1];
    xNoNaNs->size[0] = nrows;
    xNoNaNs->size[1] = 3;
    emxEnsureCapacity_real_T(xNoNaNs, irow);
    xNoNaNs_data = xNoNaNs->data;
    irow = -1;
    for (j = 0; j < n; j++) {
      if (naninfo_nNaNsInRow_data[j] == 0) {
        irow++;
        xNoNaNs_data[irow] = y_data[j];
        xNoNaNs_data[irow + xNoNaNs->size[0]] = y_data[j + b_x->size[0]];
        xNoNaNs_data[irow + xNoNaNs->size[0] * 2] =
            y_data[j + b_x->size[0] * 2];
      }
    }
    varargout_3_size =
        xzsvdc(xNoNaNs, score, latent_data, coeff_data, coeff_size);
    score_data = score->data;
    nsv = score->size[1];
    for (j = 0; j < nsv; j++) {
      irow = (nrows / 2) << 1;
      n = irow - 2;
      for (i = 0; i <= n; i += 2) {
        r = _mm_loadu_pd(&score_data[i + score->size[0] * j]);
        _mm_storeu_pd(&score_data[i + score->size[0] * j],
                      _mm_mul_pd(r, _mm_set1_pd(latent_data[j])));
      }
      for (i = irow; i < nrows; i++) {
        score_data[i + score->size[0] * j] *= latent_data[j];
      }
    }
    irow = (score->size[1] / 2) << 1;
    n = irow - 2;
    for (j = 0; j <= n; j += 2) {
      r = _mm_loadu_pd(&latent_data[0]);
      _mm_storeu_pd(&latent_data[0],
                    _mm_div_pd(_mm_mul_pd(r, r), _mm_set1_pd(DOF)));
    }
    for (j = irow; j < nsv; j++) {
      wcol = latent_data[j];
      wcol = wcol * wcol / (double)DOF;
      latent_data[j] = wcol;
    }
    if (DOF < 3) {
      nsv = score->size[1];
      if (DOF <= nsv) {
        nsv = DOF;
      }
      irow = y->size[0] * y->size[1];
      y->size[0] = nrows;
      y->size[1] = nsv;
      emxEnsureCapacity_real_T(y, irow);
      y_data = y->data;
      for (j = 0; j < nsv; j++) {
        for (i = 0; i < nrows; i++) {
          y_data[i + y->size[0] * j] = score_data[i + score->size[0] * j];
        }
      }
      varargout_3_size = nsv;
      for (j = 0; j < nsv; j++) {
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
      for (j = 0; j < irow; j++) {
        y_data[j] = score_data[j];
      }
      if (varargout_3_size - 1 >= 0) {
        memcpy(&varargout_3_data[0], &latent_data[0],
               (unsigned int)varargout_3_size * sizeof(double));
      }
      nsv = coeff_size[1];
      irow = 3 * coeff_size[1];
      if (irow - 1 >= 0) {
        memcpy(&b_coeff_data[0], &coeff_data[0],
               (unsigned int)irow * sizeof(double));
      }
    }
    n = y->size[1];
    irow = score->size[0] * score->size[1];
    score->size[0] = x->size[0];
    score->size[1] = y->size[1];
    emxEnsureCapacity_real_T(score, irow);
    score_data = score->data;
    irow = -1;
    for (i = 0; i < loop_ub; i++) {
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
    for (j = 0; j < 3; j++) {
      varargout_1_data[3 * j] = b_coeff_data[3 * j];
      irow = 3 * j + 1;
      varargout_1_data[irow] = b_coeff_data[irow];
      irow = 3 * j + 2;
      varargout_1_data[irow] = b_coeff_data[irow];
    }
    irow = varargout_2->size[0] * varargout_2->size[1];
    varargout_2->size[0] = score->size[0];
    varargout_2->size[1] = 3;
    emxEnsureCapacity_real_T(varargout_2, irow);
    xNoNaNs_data = varargout_2->data;
    for (j = 0; j < 3; j++) {
      for (i = 0; i < nrows; i++) {
        xNoNaNs_data[i + varargout_2->size[0] * j] =
            score_data[i + score->size[0] * j];
      }
    }
  } else {
    varargout_1_size[0] = 3;
    varargout_1_size[1] = nsv;
    irow = 3 * nsv;
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
    for (j = 0; j < irow; j++) {
      xNoNaNs_data[j] = score_data[j];
    }
  }
  emxFree_real_T(&score);
  nsv = varargout_1_size[1];
  for (i = 0; i < nsv; i++) {
    double absc;
    wcol = 0.0;
    xcol = 1.0;
    d = varargout_1_data[3 * i];
    absc = fabs(d);
    if (absc > 0.0) {
      wcol = absc;
      xcol = d;
      if (!rtIsNaN(d)) {
        if (d < 0.0) {
          xcol = -1.0;
        } else {
          xcol = (d > 0.0);
        }
      }
    }
    d = varargout_1_data[3 * i + 1];
    absc = fabs(d);
    if (absc > wcol) {
      wcol = absc;
      xcol = d;
      if (!rtIsNaN(d)) {
        if (d < 0.0) {
          xcol = -1.0;
        } else {
          xcol = (d > 0.0);
        }
      }
    }
    irow = 3 * i + 2;
    d = varargout_1_data[irow];
    if (fabs(d) > wcol) {
      xcol = d;
      if (!rtIsNaN(d)) {
        if (d < 0.0) {
          xcol = -1.0;
        } else {
          xcol = (d > 0.0);
        }
      }
    }
    if (xcol < 0.0) {
      __m128d r1;
      r = _mm_loadu_pd(&varargout_1_data[3 * i]);
      r1 = _mm_set1_pd(-1.0);
      _mm_storeu_pd(&varargout_1_data[3 * i], _mm_mul_pd(r, r1));
      varargout_1_data[irow] = -varargout_1_data[irow];
      irow = (nrows / 2) << 1;
      n = irow - 2;
      for (j = 0; j <= n; j += 2) {
        r = _mm_loadu_pd(&xNoNaNs_data[j + varargout_2->size[0] * i]);
        _mm_storeu_pd(&xNoNaNs_data[j + varargout_2->size[0] * i],
                      _mm_mul_pd(r, r1));
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
