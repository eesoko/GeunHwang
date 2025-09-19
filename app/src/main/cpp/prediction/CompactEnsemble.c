/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CompactEnsemble.c
 *
 * Code generation for function 'CompactEnsemble'
 *
 */

/* Include files */
#include "CompactEnsemble.h"
#include "CompactClassificationDiscriminant.h"
#include "CompactClassificationModel.h"
#include "predictOneWithCache.h"
#include "predict_exercise_index_data.h"
#include "predict_exercise_index_internal_types.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void CompactEnsemble_ensemblePredict(const double X[32], double score[6])
{
  static const double D[16] = {3.4288270160408407,
                               1.7878109925877446,
                               1.241136997066784,
                               0.8308368052781433,
                               1.2127792229395988,
                               0.13256642616028277,
                               1.1268516624651217,
                               7.3621361723689507,
                               0.0,
                               0.0,
                               1.6679238495874826,
                               0.65267069830079427,
                               1.3550610155951051,
                               28.638600306610094,
                               0.033030068921813815,
                               24.096842510754961};
  static const double S[16] = {2.3610645984360534,
                               1.2357000050350415,
                               1.1456921914664262,
                               1.0629555265421966,
                               0.947650492667453,
                               0.895890480385118,
                               0.85613165029402849,
                               0.80491892563436929,
                               0.7261898446606293,
                               0.61125615660190591,
                               0.45152159492151489,
                               0.37597619967614526,
                               0.31501308945230067,
                               0.17020118271011189,
                               0.0,
                               0.0};
  static const double b_D[16] = {
      0.8308368052781433,  0.56022048564134463,  0.42773807091983518,
      0.13256642616028277, 0.53447259241942158,  15.262218518440363,
      11.509161473433817,  11.131000098229741,   1.3550610155951051,
      28.638600306610094,  0.033030068921813815, 0.00084992667686587485,
      7.7269969729104089,  119.32937916689845,   24.096842510754961,
      0.056856379941248175};
  static const double b_S[16] = {
      1.9748139654564763,  1.8555029473646307,  1.4288940904912242,
      1.0906438299846366,  1.0435275950294078,  0.97537614595193756,
      0.93950177498072329, 0.85799715742748539, 0.68042597154860118,
      0.66401370931045844, 0.55516342997542234, 0.48098529217180708,
      0.4047175515404372,  0.27848000585351451, 0.22946865144596915,
      0.17144937168216465};
  static const double c_D[16] = {
      3.4288270160408407,   0.91209325290455545,    1.2127792229395988,
      0.56022048564134463,  0.3236507139519616,     0.33677100414213235,
      1.1268516624651217,   7.3621361723689507,     0.53447259241942158,
      11.509161473433817,   1.6679238495874826,     1.3550610155951051,
      0.033030068921813815, 0.00084992667686587485, 24.096842510754961,
      0.056856379941248175};
  static const double c_S[16] = {
      2.3255666482169577,  1.4108382353601561,  1.2217887687612285,
      1.1195790092084241,  1.0974702023290033,  0.98507685066119588,
      0.95153572312298573, 0.88515855221917461, 0.75959774599159746,
      0.71922708989141171, 0.58198395124904712, 0.4885569569162882,
      0.39618436190584033, 0.28982646159542413, 0.22430416773982995,
      0.16841554270304609};
  static const double d_D[16] = {
      3.4288270160408407,  1.7878109925877446,   1.241136997066784,
      0.3236507139519616,  1.1268516624651217,   0.53447259241942158,
      15.262218518440363,  11.509161473433817,   1.6679238495874826,
      0.65267069830079427, 0.080963196419212263, 1.3550610155951051,
      28.638600306610094,  0.033030068921813815, 0.056856379941248175,
      0.002708497303718322};
  static const double d_S[16] = {
      2.2253605717153144,  1.7885653123187786,  1.3345039994813865,
      1.1331797919400535,  1.0344941147639433,  0.96180220056579846,
      0.86904620777963559, 0.77520455671182587, 0.66716236432757559,
      0.60637236027768215, 0.52390366889862494, 0.3818982595715682,
      0.30868987025557154, 0.22289833677857396, 0.17562764830566963,
      0.15313246362179908};
  static const double e_D[16] = {1.241136997066784,
                                 1.2127792229395988,
                                 0.56022048564134463,
                                 0.42773807091983518,
                                 0.33677100414213235,
                                 0.13256642616028277,
                                 0.0,
                                 0.0,
                                 1.6679238495874826,
                                 0.080963196419212263,
                                 11.131000098229741,
                                 1.3550610155951051,
                                 0.00084992667686587485,
                                 7.7269969729104089,
                                 24.096842510754961,
                                 0.002708497303718322};
  static const double e_S[16] = {2.0602090936163422,
                                 1.5457868300380211,
                                 1.272678504395282,
                                 1.1093080402360922,
                                 0.99107807978324491,
                                 0.885271441506782,
                                 0.86703042144432829,
                                 0.74815091514303889,
                                 0.70023532481799744,
                                 0.65119317836925261,
                                 0.50051047279724037,
                                 0.41171721385472226,
                                 0.23722873778599654,
                                 0.21843190424590272,
                                 0.0,
                                 0.0};
  static const double f_D[16] = {3.4288270160408407,
                                 1.241136997066784,
                                 0.91209325290455545,
                                 0.3236507139519616,
                                 0.33677100414213235,
                                 0.13256642616028277,
                                 1.1268516624651217,
                                 7.3621361723689507,
                                 0.0,
                                 0.0,
                                 15.262218518440363,
                                 0.65267069830079427,
                                 28.638600306610094,
                                 0.033030068921813815,
                                 7.7269969729104089,
                                 0.056856379941248175};
  static const double f_S[16] = {1.8673066650166803,
                                 1.5099385009704085,
                                 1.2423473822142161,
                                 1.1478901745292467,
                                 1.091713859149011,
                                 0.99938691993655193,
                                 0.89646524599958088,
                                 0.80805316795980509,
                                 0.71983121069803968,
                                 0.6679904555485513,
                                 0.59976889067292938,
                                 0.50225287693628407,
                                 0.31840365716125446,
                                 0.21732386330204367,
                                 0.0,
                                 0.0};
  static const double g_D[16] = {3.4288270160408407,   1.7878109925877446,
                                 0.91209325290455545,  0.8308368052781433,
                                 0.42773807091983518,  0.3236507139519616,
                                 0.33677100414213235,  0.13256642616028277,
                                 0.53447259241942158,  0.0,
                                 15.262218518440363,   1.6679238495874826,
                                 0.080963196419212263, 11.131000098229741,
                                 28.638600306610094,   119.32937916689845};
  static const double g_S[16] = {2.055591395598459,   1.5751421027935968,
                                 1.342961809807673,   1.2730508135145759,
                                 1.0990957319642727,  0.88220115664644116,
                                 0.81317132390369518, 0.744182577820804,
                                 0.69177258162588062, 0.61126244779647543,
                                 0.55960312849707916, 0.46452716064700206,
                                 0.37409220758143108, 0.31811644632151465,
                                 0.2136484767637476,  0.0};
  static const double h_D[16] = {1.7878109925877446,  1.241136997066784,
                                 0.91209325290455545, 0.8308368052781433,
                                 0.42773807091983518, 0.33677100414213235,
                                 0.13256642616028277, 7.3621361723689507,
                                 0.53447259241942158, 0.0,
                                 11.509161473433817,  0.65267069830079427,
                                 1.3550610155951051,  7.7269969729104089,
                                 119.32937916689845,  24.096842510754961};
  static const double h_S[16] = {2.2946083551987013,  1.4271102725057891,
                                 1.2466935673602664,  1.1941030605182923,
                                 0.99050742633037958, 0.97113250727029488,
                                 0.88921271122139522, 0.72386586647108453,
                                 0.68848090256375138, 0.57694172479165762,
                                 0.51016103941612512, 0.42730041679573355,
                                 0.36146846558251977, 0.24407550696798849,
                                 0.19790782083332223, 0.0};
  static const double i_D[16] = {1.241136997066784,
                                 1.2127792229395988,
                                 1.1268516624651217,
                                 7.3621361723689507,
                                 0.53447259241942158,
                                 0.0,
                                 0.0,
                                 11.509161473433817,
                                 1.6679238495874826,
                                 0.65267069830079427,
                                 0.080963196419212263,
                                 1.3550610155951051,
                                 0.033030068921813815,
                                 0.00084992667686587485,
                                 119.32937916689845,
                                 24.096842510754961};
  static const double i_S[16] = {2.3879393875956172,
                                 1.3848575034413686,
                                 1.2817005556574392,
                                 0.99752357113356027,
                                 0.98542072721132623,
                                 0.93472658268775677,
                                 0.84789432852417956,
                                 0.74753843849321955,
                                 0.4761218783967669,
                                 0.45306601750112196,
                                 0.25522204924974412,
                                 0.23972165329124187,
                                 0.20349860523116775,
                                 0.15369442575919956,
                                 0.0,
                                 0.0};
  static const double j_D[16] = {3.4288270160408407,
                                 0.8308368052781433,
                                 1.2127792229395988,
                                 0.56022048564134463,
                                 0.3236507139519616,
                                 0.13256642616028277,
                                 0.0,
                                 11.509161473433817,
                                 0.080963196419212263,
                                 11.131000098229741,
                                 28.638600306610094,
                                 0.033030068921813815,
                                 119.32937916689845,
                                 24.096842510754961,
                                 0.056856379941248175,
                                 0.002708497303718322};
  static const double j_S[16] = {1.7914886314894574,  1.6339570568584851,
                                 1.4494411087525338,  1.1850492451800749,
                                 1.047973178332144,   1.0039973415180152,
                                 0.92881277617969638, 0.846139527528589,
                                 0.73132761755987785, 0.70415059802586488,
                                 0.572107891483376,   0.496237218540581,
                                 0.36754125958763983, 0.33916294136236264,
                                 0.27619344138069296, 0.0};
  static const double k_D[16] = {1.7878109925877446,     1.241136997066784,
                                 0.91209325290455545,    0.8308368052781433,
                                 0.56022048564134463,    0.42773807091983518,
                                 0.13256642616028277,    0.0,
                                 15.262218518440363,     11.509161473433817,
                                 0.65267069830079427,    0.080963196419212263,
                                 0.00084992667686587485, 7.7269969729104089,
                                 0.056856379941248175,   0.002708497303718322};
  static const double k_S[16] = {1.84258764748588,    1.5691233222214729,
                                 1.4952742201714062,  1.1355464925738421,
                                 1.0781707554000983,  0.98514021221458592,
                                 0.88995685737713082, 0.852743328567912,
                                 0.76627952036641722, 0.651398137077833,
                                 0.60332977629539075, 0.49828787084980325,
                                 0.39754209215542474, 0.34124241795558141,
                                 0.258798304448928,   0.0};
  static const char combiner[15] = {'W', 'e', 'i', 'g', 'h', 't', 'e', 'd',
                                    'A', 'v', 'e', 'r', 'a', 'g', 'e'};
  c_classreg_learning_classif_Com weak_learner;
  double cachedScore[6];
  double b_logdetsigma;
  double cachedWeights;
  double logdetsigma;
  int ii;
  int kk;
  c_CompactClassificationModel_Co(&weak_learner);
  logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma = D[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
    b_logdetsigma = S[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  for (ii = 0; ii < 6; ii++) {
    cachedScore[ii] = 0.0;
  }
  cachedWeights = predictOneWithCache(
      X, cachedScore, combiner, weak_learner.Prior,
      weak_learner.ClassLogicalIndices, weak_learner.Cost, logdetsigma, score);
  c_CompactClassificationDiscrimi(&weak_learner);
  b_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  d_CompactClassificationDiscrimi(&weak_learner);
  c_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  c_CompactClassificationModel_Co(&weak_learner);
  b_logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma += 2.0 * log(b_D[ii]);
    b_logdetsigma += 2.0 * log(b_S[ii]);
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  d_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, b_logdetsigma, score);
  c_CompactClassificationModel_Co(&weak_learner);
  b_logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma += 2.0 * log(c_D[ii]);
    b_logdetsigma += 2.0 * log(c_S[ii]);
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  e_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, b_logdetsigma, score);
  c_CompactClassificationModel_Co(&weak_learner);
  b_logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma += 2.0 * log(d_D[ii]);
    b_logdetsigma += 2.0 * log(d_S[ii]);
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  f_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, b_logdetsigma, score);
  d_CompactClassificationDiscrimi(&weak_learner);
  g_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  c_CompactClassificationModel_Co(&weak_learner);
  logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma = e_D[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
    b_logdetsigma = e_S[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  h_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, logdetsigma, score);
  e_CompactClassificationDiscrimi(&weak_learner);
  i_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  c_CompactClassificationModel_Co(&weak_learner);
  logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma = f_D[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
    b_logdetsigma = f_S[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  j_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, logdetsigma, score);
  f_CompactClassificationDiscrimi(&weak_learner);
  k_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  d_CompactClassificationDiscrimi(&weak_learner);
  l_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  c_CompactClassificationModel_Co(&weak_learner);
  logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma = g_D[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
    b_logdetsigma = g_S[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  m_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, logdetsigma, score);
  c_CompactClassificationModel_Co(&weak_learner);
  logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma = h_D[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
    b_logdetsigma = h_S[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  n_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, logdetsigma, score);
  c_CompactClassificationModel_Co(&weak_learner);
  logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma = i_D[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
    b_logdetsigma = i_S[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  o_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, logdetsigma, score);
  f_CompactClassificationDiscrimi(&weak_learner);
  p_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  c_CompactClassificationDiscrimi(&weak_learner);
  q_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  f_CompactClassificationDiscrimi(&weak_learner);
  r_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  c_CompactClassificationModel_Co(&weak_learner);
  logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma = j_D[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
    b_logdetsigma = j_S[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  s_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, logdetsigma, score);
  d_CompactClassificationDiscrimi(&weak_learner);
  t_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  d_CompactClassificationDiscrimi(&weak_learner);
  u_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  c_CompactClassificationModel_Co(&weak_learner);
  logdetsigma = 0.0;
  for (ii = 0; ii < 16; ii++) {
    b_logdetsigma = k_D[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
    b_logdetsigma = k_S[ii];
    if (b_logdetsigma > 0.0) {
      logdetsigma += 2.0 * log(b_logdetsigma);
    }
  }
  kk = 0;
  for (ii = 0; ii < 6; ii++) {
    weak_learner.Prior[ii] = 0.0;
    if (weak_learner.ClassLogicalIndices[ii]) {
      weak_learner.Prior[ii] = dv[kk];
      kk++;
    }
  }
  for (ii = 0; ii < 36; ii++) {
    weak_learner.Cost[ii] = iv[ii];
  }
  v_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, logdetsigma, score);
  g_CompactClassificationDiscrimi(&weak_learner);
  w_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  h_CompactClassificationDiscrimi(&weak_learner);
  x_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  e_CompactClassificationDiscrimi(&weak_learner);
  y_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                        weak_learner.Prior, weak_learner.ClassLogicalIndices,
                        weak_learner.Cost, weak_learner.LogDetSigma, score);
  i_CompactClassificationDiscrimi(&weak_learner);
  ab_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                         weak_learner.Prior, weak_learner.ClassLogicalIndices,
                         weak_learner.Cost, weak_learner.LogDetSigma, score);
  j_CompactClassificationDiscrimi(&weak_learner);
  bb_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                         weak_learner.Prior, weak_learner.ClassLogicalIndices,
                         weak_learner.Cost, weak_learner.LogDetSigma, score);
  d_CompactClassificationDiscrimi(&weak_learner);
  cb_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                         weak_learner.Prior, weak_learner.ClassLogicalIndices,
                         weak_learner.Cost, weak_learner.LogDetSigma, score);
  k_CompactClassificationDiscrimi(&weak_learner);
  db_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                         weak_learner.Prior, weak_learner.ClassLogicalIndices,
                         weak_learner.Cost, weak_learner.LogDetSigma, score);
  l_CompactClassificationDiscrimi(&weak_learner);
  eb_predictOneWithCache(X, cachedScore, &cachedWeights, combiner,
                         weak_learner.Prior, weak_learner.ClassLogicalIndices,
                         weak_learner.Cost, weak_learner.LogDetSigma, score);
}

/* End of code generation (CompactEnsemble.c) */
