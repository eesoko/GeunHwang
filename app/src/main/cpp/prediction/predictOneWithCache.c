/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * predictOneWithCache.c
 *
 * Code generation for function 'predictOneWithCache'
 *
 */

/* Include files */
#include "predictOneWithCache.h"
#include "CompactClassificationDiscriminant.h"
#include "CompactClassificationModel.h"
#include "ismember.h"
#include "predict_exercise_data.h"
#include "rt_nonfinite.h"
#include "updateCache.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
boolean_T
ab_predictOneWithCache(const double X[32], double cachedScore[6],
                       double *cachedWeights, const char combiner[15],
                       const double weak_learner_Prior[6],
                       const boolean_T c_weak_learner_ClassLogicalIndi[6],
                       const double weak_learner_Cost[36],
                       double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {1,  3,  6,  7,  8,  12, 14, 16,
                                       18, 20, 22, 24, 26, 28, 29, 30};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    lb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
b_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {3,  4,  6,  7,  8,  11, 12, 14,
                                       15, 16, 17, 20, 23, 26, 28, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    m_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
bb_predictOneWithCache(const double X[32], double cachedScore[6],
                       double *cachedWeights, const char combiner[15],
                       const double weak_learner_Prior[6],
                       const boolean_T c_weak_learner_ClassLogicalIndi[6],
                       const double weak_learner_Cost[36],
                       double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  3,  6,  7,  9,  12, 13, 15,
                                       17, 20, 25, 27, 28, 29, 30, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    mb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
c_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {1,  2,  3,  4,  5,  9,  14, 17,
                                       18, 19, 21, 23, 26, 27, 30, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    n_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
cb_predictOneWithCache(const double X[32], double cachedScore[6],
                       double *cachedWeights, const char combiner[15],
                       const double weak_learner_Prior[6],
                       const boolean_T c_weak_learner_ClassLogicalIndi[6],
                       const double weak_learner_Cost[36],
                       double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {1,  3,  4,  7,  9,  10, 12, 15,
                                       19, 20, 22, 23, 25, 26, 27, 30};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    nb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
d_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {4,  6,  7,  10, 13, 16, 17, 21,
                                       22, 24, 25, 26, 27, 28, 29, 30};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    o_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
db_predictOneWithCache(const double X[32], double cachedScore[6],
                       double *cachedWeights, const char combiner[15],
                       const double weak_learner_Prior[6],
                       const boolean_T c_weak_learner_ClassLogicalIndi[6],
                       const double weak_learner_Cost[36],
                       double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  1,  2,  5,  6,  11, 13, 14,
                                       17, 18, 19, 20, 26, 28, 30, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    ob_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
e_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  3,  5,  6,  8,  9,  11, 12,
                                       13, 17, 18, 22, 25, 26, 29, 30};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    p_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
eb_predictOneWithCache(const double X[32], double cachedScore[6],
                       double *cachedWeights, const char combiner[15],
                       const double weak_learner_Prior[6],
                       const boolean_T c_weak_learner_ClassLogicalIndi[6],
                       const double weak_learner_Cost[36],
                       double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  5,  6,  7,  8,  9,  10, 11,
                                       13, 15, 19, 21, 24, 25, 27, 29};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    pb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
f_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  1,  2,  8,  11, 13, 16, 17,
                                       18, 19, 20, 22, 24, 25, 30, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    q_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
g_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  3,  5,  6,  7,  8,  10, 13,
                                       19, 20, 21, 23, 24, 28, 29, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    r_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
h_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {2,  5,  6,  7,  9,  10, 14, 15,
                                       18, 20, 21, 22, 26, 27, 29, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    s_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
i_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {1,  6,  7,  8,  10, 12, 13, 16,
                                       18, 19, 20, 21, 22, 23, 25, 28};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    t_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
j_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  2,  3,  8,  9,  10, 11, 12,
                                       14, 15, 16, 19, 24, 25, 27, 30};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    u_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
k_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  1,  2,  4,  12, 14, 15, 18,
                                       21, 22, 23, 24, 25, 26, 27, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    v_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
l_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {4,  5,  6,  9,  10, 11, 12, 15,
                                       16, 17, 20, 23, 26, 28, 29, 30};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    w_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
m_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  1,  3,  4,  7,  8,  9,  10,
                                       13, 14, 16, 18, 20, 21, 24, 28};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    x_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
n_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {1,  2,  3,  4,  7,  9,  10, 12,
                                       13, 15, 17, 19, 22, 27, 28, 29};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    y_CompactClassificationDiscrimi(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
o_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {2,  5,  11, 12, 13, 14, 15, 17,
                                       18, 19, 20, 22, 25, 26, 28, 29};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    ab_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
p_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {1,  3,  5,  7,  13, 15, 16, 17,
                                       19, 21, 23, 24, 25, 27, 30, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    bb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

double predictOneWithCache(const double X[32], double cachedScore[6],
                           const char combiner[15],
                           const double weak_learner_Prior[6],
                           const boolean_T c_weak_learner_ClassLogicalIndi[6],
                           const double weak_learner_Cost[36],
                           double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  1,  2,  4,  5,  10, 11, 12,
                                       14, 15, 18, 19, 22, 24, 25, 29};
  double learnerscore[6];
  double logP_data[6];
  double mah_data[6];
  double pos[6];
  double tmp_data[6];
  double cachedWeights;
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int i1;
  int ii;
  int k;
  int trueCount;
  signed char i;
  boolean_T b[16];
  boolean_T exitg1;
  boolean_T y;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T a__1[6];
  categorical_ismember(a__1, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  y = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      y = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  trueCount = 0;
  if (!y) {
    double X_data[16];
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    qb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          cachedWeights = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          cachedWeights =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              cachedWeights;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = cachedWeights;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      cachedWeights = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (cachedWeights < d) {
          cachedWeights = d;
        }
      }
      if (rtIsInf(cachedWeights) ||
          (cachedWeights == -1.7976931348623157E+308)) {
        for (i1 = 0; i1 < 6; i1++) {
          tmp_data[logP_size_idx_0 * i1] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - cachedWeights);
        }
        cachedWeights = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          cachedWeights += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= cachedWeights;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&cachedWeights, mah_data,
                                    mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
      cachedScore[ii] = 0.0;
    }
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      mah_data[trueCount * ii] = cachedScore[ii] + learnerscore[ii];
    }
  }
  for (ii = 0; ii < 6; ii++) {
    for (i1 = 0; i1 < trueCount; i1++) {
      cachedScore[ii] = mah_data[trueCount * ii];
    }
  }
  i = 0;
  if (trueCount - 1 >= 0) {
    i = 1;
  }
  y = false;
  k = 0;
  int exitg2;
  do {
    exitg2 = 0;
    if (k < 15) {
      if (cv3[(int)combiner[k]] != cv3[(int)cv4[k]]) {
        exitg2 = 1;
      } else {
        k++;
      }
    } else {
      y = true;
      exitg2 = 1;
    }
  } while (exitg2 == 0);
  if (y) {
    if (i == 0) {
      y = true;
      k = 0;
      exitg1 = false;
      while ((!exitg1) && (k < 6)) {
        if (!(cachedScore[k] == 0.0)) {
          y = false;
          exitg1 = true;
        } else {
          k++;
        }
      }
      if (y) {
        for (ii = 0; ii < 6; ii++) {
          score[ii] = rtNaN;
        }
      } else {
        y = true;
        k = 0;
        exitg1 = false;
        while ((!exitg1) && (k < 6)) {
          if (!(cachedScore[k] < 0.0)) {
            y = false;
            exitg1 = true;
          } else {
            k++;
          }
        }
        if (y) {
          for (ii = 0; ii < 6; ii++) {
            score[ii] = rtMinusInf;
          }
        } else {
          for (ii = 0; ii < 6; ii++) {
            score[ii] = rtInf;
          }
        }
      }
    } else {
      for (ii = 0; ii < 6; ii++) {
        score[ii] = cachedScore[ii];
      }
    }
  } else {
    for (ii = 0; ii < 6; ii++) {
      score[ii] = cachedScore[ii];
    }
  }
  return i;
}

boolean_T
q_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {1,  2,  3,  5,  6,  7,  8,  9,
                                       11, 14, 21, 22, 23, 27, 29, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    cb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
r_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  4,  5,  7,  8,  9,  11, 14,
                                       15, 19, 23, 25, 28, 29, 30, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    db_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
s_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  4,  5,  6,  8,  10, 15, 17,
                                       20, 21, 24, 25, 28, 29, 30, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    eb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
t_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {4,  6,  8,  11, 12, 13, 16, 19,
                                       20, 21, 22, 23, 24, 26, 27, 28};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    fb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
u_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  2,  8,  9,  11, 12, 13, 14,
                                       16, 18, 21, 23, 24, 26, 27, 30};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    gb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
v_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {1,  2,  3,  4,  6,  7,  10, 14,
                                       16, 17, 19, 20, 26, 27, 30, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    hb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
w_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  1,  2,  10, 12, 16, 18, 20,
                                       21, 23, 24, 26, 27, 28, 29, 31};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    ib_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
x_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {0,  4,  5,  6,  7,  9,  10, 11,
                                       13, 15, 17, 18, 21, 22, 26, 29};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    jb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

boolean_T
y_predictOneWithCache(const double X[32], double cachedScore[6],
                      double *cachedWeights, const char combiner[15],
                      const double weak_learner_Prior[6],
                      const boolean_T c_weak_learner_ClassLogicalIndi[6],
                      const double weak_learner_Cost[36],
                      double weak_learner_LogDetSigma, double score[6])
{
  static const signed char b_iv[16] = {1,  2,  3,  4,  5,  6,  12, 13,
                                       14, 16, 17, 18, 22, 23, 24, 25};
  double learnerscore[6];
  double logP_data[6];
  double pos[6];
  double tmp_data[6];
  int X_size[2];
  int mah_size[2];
  int b_ii;
  int ii;
  int k;
  int trueCount;
  boolean_T b[16];
  boolean_T cached;
  boolean_T exitg1;
  boolean_T obsToUseIdx;
  for (ii = 0; ii < 6; ii++) {
    learnerscore[ii] = rtNaN;
  }
  boolean_T tf[6];
  categorical_ismember(tf, pos);
  for (ii = 0; ii < 16; ii++) {
    b[ii] = rtIsNaN(X[b_iv[ii]]);
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 16)) {
    if (b[k]) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  obsToUseIdx = !cached;
  trueCount = 0;
  if (obsToUseIdx) {
    double X_data[16];
    double mah_data[6];
    double sumP_data;
    int logP_size_idx_0;
    trueCount = 1;
    mah_size[0] = 1;
    mah_size[1] = 6;
    for (ii = 0; ii < 6; ii++) {
      mah_data[ii] = rtNaN;
    }
    X_size[0] = 1;
    X_size[1] = 16;
    for (ii = 0; ii < 16; ii++) {
      X_data[ii] = X[b_iv[ii]];
    }
    kb_CompactClassificationDiscrim(
        mah_data, mah_size, c_weak_learner_ClassLogicalIndi, X_data, X_size);
    logP_size_idx_0 = (signed char)mah_size[0];
    k = (signed char)mah_size[0] * 6;
    if (k - 1 >= 0) {
      memset(&logP_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (ii = 0; ii < 6; ii++) {
      if (weak_learner_Prior[ii] == 0.0) {
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = rtMinusInf;
        }
      } else {
        int loop_ub;
        loop_ub = mah_size[0];
        for (b_ii = 0; b_ii < loop_ub; b_ii++) {
          sumP_data = mah_data[mah_size[0] * ii] / 2.0;
        }
        if (mah_size[0] != 0) {
          sumP_data =
              (log(weak_learner_Prior[ii]) - weak_learner_LogDetSigma / 2.0) -
              sumP_data;
        }
        if (logP_size_idx_0 - 1 >= 0) {
          logP_data[logP_size_idx_0 * ii] = sumP_data;
        }
      }
    }
    X_size[0] = (signed char)mah_size[0];
    X_size[1] = 6;
    if (k - 1 >= 0) {
      memset(&tmp_data[0], 0, (unsigned int)k * sizeof(double));
    }
    for (b_ii = 0; b_ii < logP_size_idx_0; b_ii++) {
      sumP_data = logP_data[0];
      for (ii = 0; ii < 6; ii++) {
        double d;
        d = logP_data[logP_size_idx_0 * ii];
        if (sumP_data < d) {
          sumP_data = d;
        }
      }
      if (rtIsInf(sumP_data) || (sumP_data == -1.7976931348623157E+308)) {
        for (ii = 0; ii < 6; ii++) {
          tmp_data[logP_size_idx_0 * ii] = rtNaN;
        }
      } else {
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] = exp(logP_data[k] - sumP_data);
        }
        sumP_data = tmp_data[0];
        for (ii = 0; ii < 5; ii++) {
          sumP_data += tmp_data[logP_size_idx_0 * (ii + 1)];
        }
        for (ii = 0; ii < 6; ii++) {
          k = logP_size_idx_0 * ii;
          tmp_data[k] /= sumP_data;
        }
      }
    }
    c_CompactClassificationModel_mi(weak_learner_Cost, tmp_data, X_size,
                                    (double *)&sumP_data, mah_data, mah_size);
  } else {
    X_size[0] = 0;
  }
  for (ii = 0; ii < 6; ii++) {
    if (trueCount - 1 >= 0) {
      learnerscore[(int)pos[ii] - 1] = tmp_data[X_size[0] * ii];
    }
  }
  cached = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (rtIsNaN(cachedScore[k])) {
      cached = true;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (cached) {
    for (ii = 0; ii < 6; ii++) {
      if (rtIsNaN(cachedScore[ii]) && obsToUseIdx) {
        cachedScore[ii] = 0.0;
      }
    }
  }
  cached = false;
  updateCache(learnerscore, cachedScore, cachedWeights, &cached, combiner,
              obsToUseIdx, score);
  return true;
}

/* End of code generation (predictOneWithCache.c) */
