//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// pca.cpp
//
// Code generation for function 'pca'
//

// Include files
#include "pca.h"
#include "rt_nonfinite.h"
#include "xzsvdc.h"
#include "coder_array.h"
#include <algorithm>
#include <cmath>
#include <cstring>

// Function Definitions
namespace coder {
int pca(const array<double, 2U> &x, double varargout_1_data[],
        int varargout_1_size[2], array<double, 2U> &varargout_2,
        double varargout_3_data[])
{
  array<double, 2U> b_x;
  array<double, 2U> score;
  array<double, 2U> xNoNaNs;
  array<double, 2U> y;
  array<int, 1U> naninfo_nNaNsInRow;
  array<bool, 2U> naninfo_isNaN;
  double b_coeff_data[9];
  double coeff_data[9];
  double latent_data[3];
  double mu[3];
  double sgn;
  double wcol;
  double xcol;
  int coeff_size[2];
  int DOF;
  int b_nsv;
  int loop_ub;
  int minDOFnsv;
  int ns1;
  int nsv;
  int p;
  int varargout_3_size;
  bool noNaNs;
  loop_ub = x.size(0);
  b_x.set_size(x.size(0), 3);
  nsv = x.size(0) * 3;
  for (int j{0}; j < nsv; j++) {
    b_x[j] = x[j];
  }
  b_nsv = 0;
  p = 0;
  naninfo_nNaNsInRow.set_size(x.size(0));
  if (loop_ub - 1 >= 0) {
    std::memset(&naninfo_nNaNsInRow[0], 0,
                static_cast<unsigned int>(loop_ub) * sizeof(int));
  }
  naninfo_isNaN.set_size(x.size(0), 3);
  for (int j{0}; j < nsv; j++) {
    naninfo_isNaN[j] = std::isnan(x[j]);
  }
  for (int j{0}; j < 3; j++) {
    for (int i{0}; i < loop_ub; i++) {
      if (naninfo_isNaN[i + naninfo_isNaN.size(0) * j]) {
        naninfo_nNaNsInRow[i] = naninfo_nNaNsInRow[i] + 1;
        b_nsv++;
      }
    }
  }
  for (int j{0}; j < loop_ub; j++) {
    if (naninfo_nNaNsInRow[j] > 0) {
      p++;
    }
  }
  noNaNs = (b_nsv <= 0);
  ns1 = x.size(0) - p;
  DOF = ns1;
  if (ns1 >= 1) {
    DOF = ns1 - 1;
  }
  if (!noNaNs) {
    for (int j{0}; j < 3; j++) {
      wcol = 0.0;
      xcol = 0.0;
      for (int i{0}; i < loop_ub; i++) {
        sgn = x[i + x.size(0) * j];
        if (!std::isnan(sgn)) {
          wcol++;
          xcol += sgn;
        }
      }
      mu[j] = xcol / wcol;
    }
  } else {
    for (int j{0}; j < 3; j++) {
      wcol = 0.0;
      xcol = 0.0;
      for (int i{0}; i < loop_ub; i++) {
        wcol++;
        xcol += x[i + x.size(0) * j];
      }
      mu[j] = xcol / wcol;
    }
  }
  for (int j{0}; j < 3; j++) {
    for (int i{0}; i < loop_ub; i++) {
      nsv = i + b_x.size(0) * j;
      b_x[nsv] = b_x[nsv] - mu[j];
    }
  }
  if (noNaNs) {
    varargout_3_size = internal::reflapack::xzsvdc(b_x, y, latent_data,
                                                   coeff_data, coeff_size);
    nsv = y.size(1);
    for (int j{0}; j < nsv; j++) {
      for (int i{0}; i < loop_ub; i++) {
        b_nsv = i + y.size(0) * j;
        y[b_nsv] = y[b_nsv] * latent_data[j];
      }
      wcol = latent_data[j];
      wcol = wcol * wcol / static_cast<double>(DOF);
      latent_data[j] = wcol;
    }
    if (DOF < 3) {
      minDOFnsv = y.size(1);
      if (DOF <= minDOFnsv) {
        minDOFnsv = DOF;
      }
      score.set_size(x.size(0), minDOFnsv);
      varargout_3_size = minDOFnsv;
      for (int i{0}; i < minDOFnsv; i++) {
        for (int j{0}; j < loop_ub; j++) {
          score[j + score.size(0) * i] = y[j + y.size(0) * i];
        }
        varargout_3_data[i] = latent_data[i];
        b_coeff_data[3 * i] = coeff_data[3 * i];
        nsv = 3 * i + 1;
        b_coeff_data[nsv] = coeff_data[nsv];
        nsv = 3 * i + 2;
        b_coeff_data[nsv] = coeff_data[nsv];
      }
    } else {
      score.set_size(y.size(0), y.size(1));
      nsv = y.size(0) * y.size(1);
      if (nsv - 1 >= 0) {
        std::copy(&y[0], &y[nsv], &score[0]);
      }
      if (varargout_3_size - 1 >= 0) {
        std::copy(&latent_data[0], &latent_data[varargout_3_size],
                  &varargout_3_data[0]);
      }
      minDOFnsv = coeff_size[1];
      nsv = 3 * coeff_size[1];
      if (nsv - 1 >= 0) {
        std::copy(&coeff_data[0], &coeff_data[nsv], &b_coeff_data[0]);
      }
    }
  } else {
    xNoNaNs.set_size(ns1, 3);
    nsv = -1;
    for (int j{0}; j < loop_ub; j++) {
      if (naninfo_nNaNsInRow[j] == 0) {
        nsv++;
        xNoNaNs[nsv] = b_x[j];
        xNoNaNs[nsv + xNoNaNs.size(0)] = b_x[j + b_x.size(0)];
        xNoNaNs[nsv + xNoNaNs.size(0) * 2] = b_x[j + b_x.size(0) * 2];
      }
    }
    varargout_3_size = internal::reflapack::xzsvdc(xNoNaNs, score, latent_data,
                                                   coeff_data, coeff_size);
    b_nsv = score.size(1);
    for (int j{0}; j < b_nsv; j++) {
      for (int i{0}; i < ns1; i++) {
        p = i + score.size(0) * j;
        score[p] = score[p] * latent_data[j];
      }
      wcol = latent_data[j];
      wcol = wcol * wcol / static_cast<double>(DOF);
      latent_data[j] = wcol;
    }
    if (DOF < 3) {
      minDOFnsv = score.size(1);
      if (DOF <= minDOFnsv) {
        minDOFnsv = DOF;
      }
      y.set_size(ns1, minDOFnsv);
      varargout_3_size = minDOFnsv;
      for (int i{0}; i < minDOFnsv; i++) {
        for (int j{0}; j < ns1; j++) {
          y[j + y.size(0) * i] = score[j + score.size(0) * i];
        }
        varargout_3_data[i] = latent_data[i];
        b_coeff_data[3 * i] = coeff_data[3 * i];
        b_nsv = 3 * i + 1;
        b_coeff_data[b_nsv] = coeff_data[b_nsv];
        b_nsv = 3 * i + 2;
        b_coeff_data[b_nsv] = coeff_data[b_nsv];
      }
    } else {
      y.set_size(score.size(0), score.size(1));
      nsv = score.size(0) * score.size(1);
      if (nsv - 1 >= 0) {
        std::copy(&score[0], &score[nsv], &y[0]);
      }
      if (varargout_3_size - 1 >= 0) {
        std::copy(&latent_data[0], &latent_data[varargout_3_size],
                  &varargout_3_data[0]);
      }
      minDOFnsv = coeff_size[1];
      nsv = 3 * coeff_size[1];
      if (nsv - 1 >= 0) {
        std::copy(&coeff_data[0], &coeff_data[nsv], &b_coeff_data[0]);
      }
    }
    p = y.size(1);
    score.set_size(x.size(0), y.size(1));
    nsv = -1;
    for (int i{0}; i < loop_ub; i++) {
      if (naninfo_nNaNsInRow[i] > 0) {
        for (int j{0}; j < p; j++) {
          score[i + score.size(0) * j] = rtNaN;
        }
      } else {
        nsv++;
        for (int j{0}; j < p; j++) {
          score[i + score.size(0) * j] = y[nsv + y.size(0) * j];
        }
      }
    }
  }
  ns1 = score.size(0);
  if (DOF > 3) {
    varargout_1_size[0] = 3;
    varargout_1_size[1] = 3;
    for (int j{0}; j < 3; j++) {
      varargout_1_data[3 * j] = b_coeff_data[3 * j];
      nsv = 3 * j + 1;
      varargout_1_data[nsv] = b_coeff_data[nsv];
      nsv = 3 * j + 2;
      varargout_1_data[nsv] = b_coeff_data[nsv];
    }
    varargout_2.set_size(score.size(0), 3);
    for (int j{0}; j < 3; j++) {
      for (int i{0}; i < ns1; i++) {
        varargout_2[i + varargout_2.size(0) * j] = score[i + score.size(0) * j];
      }
    }
  } else {
    varargout_1_size[0] = 3;
    varargout_1_size[1] = minDOFnsv;
    nsv = 3 * minDOFnsv;
    if (nsv - 1 >= 0) {
      std::copy(&b_coeff_data[0], &b_coeff_data[nsv], &varargout_1_data[0]);
    }
    varargout_2.set_size(score.size(0), score.size(1));
    nsv = score.size(0) * score.size(1);
    for (int j{0}; j < nsv; j++) {
      varargout_2[j] = score[j];
    }
  }
  p = varargout_1_size[1];
  for (int i{0}; i < p; i++) {
    double d;
    double d1;
    xcol = 0.0;
    sgn = 1.0;
    d = varargout_1_data[3 * i];
    wcol = std::abs(d);
    if (wcol > 0.0) {
      xcol = wcol;
      sgn = d;
      if (!std::isnan(d)) {
        if (d < 0.0) {
          sgn = -1.0;
        } else {
          sgn = (d > 0.0);
        }
      }
    }
    nsv = 3 * i + 1;
    d1 = varargout_1_data[nsv];
    wcol = std::abs(d1);
    if (wcol > xcol) {
      xcol = wcol;
      sgn = d1;
      if (!std::isnan(d1)) {
        if (d1 < 0.0) {
          sgn = -1.0;
        } else {
          sgn = (d1 > 0.0);
        }
      }
    }
    b_nsv = 3 * i + 2;
    wcol = varargout_1_data[b_nsv];
    if (std::abs(wcol) > xcol) {
      sgn = wcol;
      if (!std::isnan(wcol)) {
        if (wcol < 0.0) {
          sgn = -1.0;
        } else {
          sgn = (wcol > 0.0);
        }
      }
    }
    if (sgn < 0.0) {
      d = -d;
      varargout_1_data[3 * i] = d;
      d1 = -d1;
      varargout_1_data[nsv] = d1;
      wcol = -wcol;
      varargout_1_data[b_nsv] = wcol;
      for (int j{0}; j < ns1; j++) {
        nsv = j + varargout_2.size(0) * i;
        varargout_2[nsv] = -varargout_2[nsv];
      }
    }
  }
  return varargout_3_size;
}

} // namespace coder

// End of code generation (pca.cpp)
