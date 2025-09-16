/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlascl.c
 *
 * Code generation for function 'xzlascl'
 *
 */

/* Include files */
#include "xzlascl.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>

/* Function Definitions */
void b_xzlascl(double cfrom, double cto, int m, double A_data[], int lda)
{
  double cfromc;
  double ctoc;
  int b_i;
  int j;
  boolean_T notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    double cfrom1;
    double cto1;
    double mul;
    int i;
    int scalarLB;
    int vectorUB;
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (fabs(cto1) > fabs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }
    i = (unsigned char)m;
    scalarLB = ((unsigned char)m >> 1) << 1;
    vectorUB = scalarLB - 2;
    for (j = 0; j < 3; j++) {
      int i1;
      int offset;
      offset = j * lda - 1;
      for (b_i = 0; b_i <= vectorUB; b_i += 2) {
        __m128d r;
        i1 = (offset + b_i) + 1;
        r = _mm_loadu_pd(&A_data[i1]);
        _mm_storeu_pd(&A_data[i1], _mm_mul_pd(r, _mm_set1_pd(mul)));
      }
      for (b_i = scalarLB; b_i < i; b_i++) {
        i1 = (offset + b_i) + 1;
        A_data[i1] *= mul;
      }
    }
  }
}

void xzlascl(double cfrom, double cto, int m, double A_data[])
{
  double cfromc;
  double ctoc;
  int b_i;
  boolean_T notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    double cfrom1;
    double cto1;
    double mul;
    int i;
    int scalarLB;
    int vectorUB;
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (fabs(cto1) > fabs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }
    i = (unsigned char)m;
    scalarLB = ((unsigned char)m >> 1) << 1;
    vectorUB = scalarLB - 2;
    for (b_i = 0; b_i <= vectorUB; b_i += 2) {
      __m128d r;
      r = _mm_loadu_pd(&A_data[b_i]);
      _mm_storeu_pd(&A_data[b_i], _mm_mul_pd(r, _mm_set1_pd(mul)));
    }
    for (b_i = scalarLB; b_i < i; b_i++) {
      A_data[b_i] *= mul;
    }
  }
}

/* End of code generation (xzlascl.c) */
