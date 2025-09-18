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
#include "rt_nonfinite.h"
#include "xzsvdc.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Definitions */
int pca(const double x_data[], const int x_size[2], double varargout_1_data[],
        int varargout_1_size[2], double varargout_2_data[],
        int varargout_2_size[2], double varargout_3_data[])
{
  __m128d r;
  double b_x_data[450];
  double c_x_data[450];
  double score_data[450];
  double y_data[450];
  double b_coeff_data[9];
  double coeff_data[9];
  double mu[3];
  double d;
  double wcol;
  double xcol;
  int naninfo_nNaNsInRow_data[150];
  int b_x_size[2];
  int y_size[2];
  int DOF;
  int DOF_tmp;
  int b_irow;
  int i;
  int irow;
  int j;
  int loop_ub;
  int minDOFnsv;
  int ns1;
  int nsv;
  int varargout_3_size;
  boolean_T naninfo_isNaN_data[450];
  boolean_T noNaNs;
  loop_ub = x_size[0];
  ns1 = x_size[0];
  irow = x_size[0] * 3;
  if (irow - 1 >= 0) {
    memcpy(&b_x_data[0], &x_data[0], (unsigned int)irow * sizeof(double));
  }
  nsv = 0;
  minDOFnsv = 0;
  if (loop_ub - 1 >= 0) {
    memset(&naninfo_nNaNsInRow_data[0], 0, (unsigned int)loop_ub * sizeof(int));
  }
  b_irow = x_size[0];
  for (j = 0; j < irow; j++) {
    naninfo_isNaN_data[j] = rtIsNaN(x_data[j]);
  }
  for (j = 0; j < 3; j++) {
    for (i = 0; i < loop_ub; i++) {
      if (naninfo_isNaN_data[i + b_irow * j]) {
        naninfo_nNaNsInRow_data[i]++;
        nsv++;
      }
    }
  }
  for (j = 0; j < loop_ub; j++) {
    if (naninfo_nNaNsInRow_data[j] > 0) {
      minDOFnsv++;
    }
  }
  noNaNs = (nsv <= 0);
  DOF_tmp = x_size[0] - minDOFnsv;
  DOF = DOF_tmp;
  if (DOF_tmp >= 1) {
    DOF = DOF_tmp - 1;
  }
  if (!noNaNs) {
    for (j = 0; j < 3; j++) {
      wcol = 0.0;
      xcol = 0.0;
      for (i = 0; i < loop_ub; i++) {
        d = x_data[i + x_size[0] * j];
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
        xcol += x_data[i + x_size[0] * j];
      }
      mu[j] = xcol / wcol;
    }
  }
  minDOFnsv = (x_size[0] / 2) << 1;
  b_irow = minDOFnsv - 2;
  for (j = 0; j < 3; j++) {
    for (i = 0; i <= b_irow; i += 2) {
      irow = i + ns1 * j;
      r = _mm_loadu_pd(&b_x_data[irow]);
      _mm_storeu_pd(&b_x_data[irow], _mm_sub_pd(r, _mm_set1_pd(mu[j])));
    }
    for (i = minDOFnsv; i < loop_ub; i++) {
      irow = i + ns1 * j;
      b_x_data[irow] -= mu[j];
    }
  }
  if (noNaNs) {
    int coeff_size[2];
    b_x_size[0] = x_size[0];
    b_x_size[1] = 3;
    irow = x_size[0] * 3 - 1;
    if (irow >= 0) {
      memcpy(&c_x_data[0], &b_x_data[0],
             (unsigned int)(irow + 1) * sizeof(double));
    }
    varargout_3_size =
        xzsvdc(c_x_data, b_x_size, y_data, y_size, mu, coeff_data, coeff_size);
    nsv = y_size[1];
    for (j = 0; j < nsv; j++) {
      irow = minDOFnsv - 2;
      for (i = 0; i <= irow; i += 2) {
        b_irow = i + y_size[0] * j;
        r = _mm_loadu_pd(&y_data[b_irow]);
        _mm_storeu_pd(&y_data[b_irow], _mm_mul_pd(r, _mm_set1_pd(mu[j])));
      }
      for (i = minDOFnsv; i < ns1; i++) {
        irow = i + y_size[0] * j;
        y_data[irow] *= mu[j];
      }
    }
    irow = (y_size[1] / 2) << 1;
    b_irow = irow - 2;
    for (j = 0; j <= b_irow; j += 2) {
      r = _mm_loadu_pd(&mu[0]);
      _mm_storeu_pd(&mu[0], _mm_div_pd(_mm_mul_pd(r, r), _mm_set1_pd(DOF)));
    }
    for (j = irow; j < nsv; j++) {
      wcol = mu[j];
      wcol = wcol * wcol / (double)DOF;
      mu[j] = wcol;
    }
    if (DOF < 3) {
      minDOFnsv = y_size[1];
      if (DOF <= minDOFnsv) {
        minDOFnsv = DOF;
      }
      varargout_2_size[0] = x_size[0];
      varargout_2_size[1] = minDOFnsv;
      varargout_3_size = minDOFnsv;
      for (i = 0; i < minDOFnsv; i++) {
        for (j = 0; j < ns1; j++) {
          score_data[j + ns1 * i] = y_data[j + y_size[0] * i];
        }
        varargout_3_data[i] = mu[i];
        b_coeff_data[3 * i] = coeff_data[3 * i];
        irow = 3 * i + 1;
        b_coeff_data[irow] = coeff_data[irow];
        irow = 3 * i + 2;
        b_coeff_data[irow] = coeff_data[irow];
      }
    } else {
      varargout_2_size[0] = y_size[0];
      varargout_2_size[1] = y_size[1];
      irow = y_size[0] * y_size[1];
      if (irow - 1 >= 0) {
        memcpy(&score_data[0], &y_data[0], (unsigned int)irow * sizeof(double));
      }
      if (varargout_3_size - 1 >= 0) {
        memcpy(&varargout_3_data[0], &mu[0],
               (unsigned int)varargout_3_size * sizeof(double));
      }
      minDOFnsv = coeff_size[1];
      irow = 3 * coeff_size[1];
      if (irow - 1 >= 0) {
        memcpy(&b_coeff_data[0], &coeff_data[0],
               (unsigned int)irow * sizeof(double));
      }
    }
  } else {
    int coeff_size[2];
    b_x_size[0] = DOF_tmp;
    b_x_size[1] = 3;
    irow = -1;
    for (j = 0; j < ns1; j++) {
      if (naninfo_nNaNsInRow_data[j] == 0) {
        irow++;
        c_x_data[irow] = b_x_data[j];
        c_x_data[irow + DOF_tmp] = b_x_data[j + ns1];
        c_x_data[irow + DOF_tmp * 2] = b_x_data[j + ns1 * 2];
      }
    }
    varargout_3_size = xzsvdc(c_x_data, b_x_size, score_data, varargout_2_size,
                              mu, coeff_data, coeff_size);
    ns1 = varargout_2_size[1];
    for (j = 0; j < ns1; j++) {
      minDOFnsv = (DOF_tmp / 2) << 1;
      irow = minDOFnsv - 2;
      for (i = 0; i <= irow; i += 2) {
        nsv = i + varargout_2_size[0] * j;
        r = _mm_loadu_pd(&score_data[nsv]);
        _mm_storeu_pd(&score_data[nsv], _mm_mul_pd(r, _mm_set1_pd(mu[j])));
      }
      for (i = minDOFnsv; i < DOF_tmp; i++) {
        irow = i + varargout_2_size[0] * j;
        score_data[irow] *= mu[j];
      }
    }
    irow = (varargout_2_size[1] / 2) << 1;
    b_irow = irow - 2;
    for (j = 0; j <= b_irow; j += 2) {
      r = _mm_loadu_pd(&mu[0]);
      _mm_storeu_pd(&mu[0], _mm_div_pd(_mm_mul_pd(r, r), _mm_set1_pd(DOF)));
    }
    for (j = irow; j < ns1; j++) {
      wcol = mu[j];
      wcol = wcol * wcol / (double)DOF;
      mu[j] = wcol;
    }
    if (DOF < 3) {
      minDOFnsv = varargout_2_size[1];
      if (DOF <= minDOFnsv) {
        minDOFnsv = DOF;
      }
      y_size[0] = DOF_tmp;
      y_size[1] = minDOFnsv;
      varargout_3_size = minDOFnsv;
      for (i = 0; i < minDOFnsv; i++) {
        for (j = 0; j < DOF_tmp; j++) {
          y_data[j + DOF_tmp * i] = score_data[j + varargout_2_size[0] * i];
        }
        varargout_3_data[i] = mu[i];
        b_coeff_data[3 * i] = coeff_data[3 * i];
        irow = 3 * i + 1;
        b_coeff_data[irow] = coeff_data[irow];
        irow = 3 * i + 2;
        b_coeff_data[irow] = coeff_data[irow];
      }
    } else {
      y_size[0] = varargout_2_size[0];
      y_size[1] = varargout_2_size[1];
      irow = varargout_2_size[0] * varargout_2_size[1];
      if (irow - 1 >= 0) {
        memcpy(&y_data[0], &score_data[0], (unsigned int)irow * sizeof(double));
      }
      if (varargout_3_size - 1 >= 0) {
        memcpy(&varargout_3_data[0], &mu[0],
               (unsigned int)varargout_3_size * sizeof(double));
      }
      minDOFnsv = coeff_size[1];
      irow = 3 * coeff_size[1];
      if (irow - 1 >= 0) {
        memcpy(&b_coeff_data[0], &coeff_data[0],
               (unsigned int)irow * sizeof(double));
      }
    }
    irow = y_size[1];
    varargout_2_size[0] = x_size[0];
    varargout_2_size[1] = y_size[1];
    b_irow = -1;
    for (i = 0; i < loop_ub; i++) {
      if (naninfo_nNaNsInRow_data[i] > 0) {
        for (j = 0; j < irow; j++) {
          score_data[i + varargout_2_size[0] * j] = rtNaN;
        }
      } else {
        b_irow++;
        for (j = 0; j < irow; j++) {
          score_data[i + varargout_2_size[0] * j] =
              y_data[b_irow + y_size[0] * j];
        }
      }
    }
  }
  ns1 = varargout_2_size[0];
  if (DOF > 3) {
    varargout_1_size[0] = 3;
    varargout_1_size[1] = 3;
    varargout_2_size[1] = 3;
    for (i = 0; i < 3; i++) {
      varargout_1_data[3 * i] = b_coeff_data[3 * i];
      irow = 3 * i + 1;
      varargout_1_data[irow] = b_coeff_data[irow];
      irow = 3 * i + 2;
      varargout_1_data[irow] = b_coeff_data[irow];
      for (j = 0; j < ns1; j++) {
        irow = varargout_2_size[0] * i;
        varargout_2_data[j + irow] = score_data[j + irow];
      }
    }
  } else {
    varargout_1_size[0] = 3;
    varargout_1_size[1] = minDOFnsv;
    irow = 3 * minDOFnsv;
    if (irow - 1 >= 0) {
      memcpy(&varargout_1_data[0], &b_coeff_data[0],
             (unsigned int)irow * sizeof(double));
    }
    irow = varargout_2_size[0] * varargout_2_size[1];
    if (irow - 1 >= 0) {
      memcpy(&varargout_2_data[0], &score_data[0],
             (unsigned int)irow * sizeof(double));
    }
  }
  minDOFnsv = varargout_1_size[1];
  for (i = 0; i < minDOFnsv; i++) {
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
      nsv = (ns1 / 2) << 1;
      irow = nsv - 2;
      for (j = 0; j <= irow; j += 2) {
        b_irow = j + varargout_2_size[0] * i;
        r = _mm_loadu_pd(&varargout_2_data[b_irow]);
        _mm_storeu_pd(&varargout_2_data[b_irow], _mm_mul_pd(r, r1));
      }
      for (j = nsv; j < ns1; j++) {
        irow = j + varargout_2_size[0] * i;
        varargout_2_data[irow] = -varargout_2_data[irow];
      }
    }
  }
  return varargout_3_size;
}

/* End of code generation (pca.c) */
