/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzsvdc.c
 *
 * Code generation for function 'xzsvdc'
 *
 */

/* Include files */
#include "xzsvdc.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "xaxpy.h"
#include "xdotc.h"
#include "xnrm2.h"
#include "xrot.h"
#include "xrotg.h"
#include "xswap.h"
#include "xzlangeM.h"
#include "xzlascl.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Definitions */
int xzsvdc(emxArray_real_T *A, emxArray_real_T *U, double S_data[],
           double V_data[], int V_size[2])
{
  emxArray_real_T *work;
  double Vf[9];
  double e[3];
  double s_data[3];
  double anrm;
  double b;
  double cscale;
  double f;
  double nrm;
  double sm;
  double sqds;
  double *A_data;
  double *U_data;
  double *work_data;
  int S_size;
  int j;
  int k;
  int n;
  int ns;
  int q;
  boolean_T doscale;
  A_data = A->data;
  n = A->size[0];
  ns = A->size[0] + 1;
  if (ns > 3) {
    ns = 3;
  }
  S_size = A->size[0];
  if (S_size > 3) {
    S_size = 3;
  }
  if (ns - 1 >= 0) {
    memset(&s_data[0], 0, (unsigned int)ns * sizeof(double));
  }
  e[0] = 0.0;
  e[1] = 0.0;
  e[2] = 0.0;
  emxInit_real_T(&work, 1);
  ns = work->size[0];
  work->size[0] = A->size[0];
  emxEnsureCapacity_real_T(work, ns);
  work_data = work->data;
  for (j = 0; j < n; j++) {
    work_data[j] = 0.0;
  }
  ns = U->size[0] * U->size[1];
  U->size[0] = A->size[0];
  U->size[1] = S_size;
  emxEnsureCapacity_real_T(U, ns);
  U_data = U->data;
  ns = A->size[0] * S_size;
  for (j = 0; j < ns; j++) {
    U_data[j] = 0.0;
  }
  memset(&Vf[0], 0, 9U * sizeof(double));
  doscale = false;
  cscale = 0.0;
  anrm = xzlangeM(A);
  if (A->size[0] == 0) {
    Vf[0] = 1.0;
    Vf[4] = 1.0;
    Vf[8] = 1.0;
  } else {
    __m128d r;
    double rt;
    double snorm;
    int iter;
    int ix;
    int iy;
    int m;
    int nct;
    int nctp1;
    int nmq;
    int nrt;
    int qp1;
    int qq;
    boolean_T guard1;
    cscale = anrm;
    guard1 = false;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      doscale = true;
      cscale = 6.7178761075670888E-139;
      guard1 = true;
    } else if (anrm > 1.4885657073574029E+138) {
      doscale = true;
      cscale = 1.4885657073574029E+138;
      guard1 = true;
    }
    if (guard1) {
      b_xzlascl(anrm, cscale, A->size[0], A, A->size[0]);
      A_data = A->data;
    }
    nrt = (A->size[0] >= 1);
    if (A->size[0] >= 1) {
      nct = A->size[0] - 1;
    } else {
      nct = 0;
    }
    if (nct > 3) {
      nct = 3;
    }
    nctp1 = nct + 1;
    if (nct >= nrt) {
      iter = nct;
    } else {
      iter = 1;
    }
    for (q = 0; q < iter; q++) {
      boolean_T apply_transform;
      qp1 = q + 2;
      qq = (q + n * q) + 1;
      nmq = (n - q) - 1;
      apply_transform = false;
      if (q + 1 <= nct) {
        nrm = xnrm2(nmq + 1, A, qq);
        if (nrm > 0.0) {
          apply_transform = true;
          if (A_data[qq - 1] < 0.0) {
            nrm = -nrm;
          }
          s_data[q] = nrm;
          if (fabs(nrm) >= 1.0020841800044864E-292) {
            nrm = 1.0 / nrm;
            ns = qq + nmq;
            ix = ((((ns - qq) + 1) / 2) << 1) + qq;
            iy = ix - 2;
            for (j = qq; j <= iy; j += 2) {
              r = _mm_loadu_pd(&A_data[j - 1]);
              _mm_storeu_pd(&A_data[j - 1], _mm_mul_pd(_mm_set1_pd(nrm), r));
            }
            for (j = ix; j <= ns; j++) {
              A_data[j - 1] *= nrm;
            }
          } else {
            ns = qq + nmq;
            ix = ((((ns - qq) + 1) / 2) << 1) + qq;
            iy = ix - 2;
            for (j = qq; j <= iy; j += 2) {
              r = _mm_loadu_pd(&A_data[j - 1]);
              _mm_storeu_pd(&A_data[j - 1],
                            _mm_div_pd(r, _mm_set1_pd(s_data[q])));
            }
            for (j = ix; j <= ns; j++) {
              A_data[j - 1] /= s_data[q];
            }
          }
          A_data[qq - 1]++;
          s_data[q] = -s_data[q];
        } else {
          s_data[q] = 0.0;
        }
      }
      for (j = qp1; j < 4; j++) {
        ns = q + n * (j - 1);
        if (apply_transform) {
          nrm = 0.0;
          if (nmq >= 0) {
            for (k = 0; k <= nmq; k++) {
              nrm += A_data[(qq + k) - 1] * A_data[ns + k];
            }
          }
          nrm = -(nrm / A_data[q + A->size[0] * q]);
          xaxpy(nmq + 1, nrm, qq, A, ns + 1);
          A_data = A->data;
        }
        e[j - 1] = A_data[ns];
      }
      if (q + 1 <= nct) {
        for (j = q + 1; j <= n; j++) {
          U_data[(j + U->size[0] * q) - 1] = A_data[(j + A->size[0] * q) - 1];
        }
      }
      if (q + 1 <= nrt) {
        nrm = b_xnrm2(e);
        if (nrm == 0.0) {
          e[0] = 0.0;
        } else {
          if (e[1] < 0.0) {
            e[0] = -nrm;
          } else {
            e[0] = nrm;
          }
          nrm = e[0];
          if (fabs(e[0]) >= 1.0020841800044864E-292) {
            nrm = 1.0 / e[0];
            ns = ((((2 - q) / 2) << 1) + q) + 2;
            iy = ns - 2;
            for (j = qp1; j <= iy; j += 2) {
              r = _mm_loadu_pd(&e[j - 1]);
              _mm_storeu_pd(&e[j - 1], _mm_mul_pd(_mm_set1_pd(nrm), r));
            }
            for (j = ns; j < 4; j++) {
              e[j - 1] *= nrm;
            }
          } else {
            ns = ((((2 - q) / 2) << 1) + q) + 2;
            iy = ns - 2;
            for (j = qp1; j <= iy; j += 2) {
              r = _mm_loadu_pd(&e[j - 1]);
              _mm_storeu_pd(&e[j - 1], _mm_div_pd(r, _mm_set1_pd(nrm)));
            }
            for (j = ns; j < 4; j++) {
              e[j - 1] /= nrm;
            }
          }
          e[1]++;
          e[0] = -e[0];
          if (n >= 2) {
            for (j = qp1; j <= n; j++) {
              work_data[j - 1] = 0.0;
            }
            for (j = qp1; j < 4; j++) {
              c_xaxpy(nmq, e[j - 1], A, n * (j - 1) + 2, work);
              work_data = work->data;
            }
            for (j = qp1; j < 4; j++) {
              d_xaxpy(nmq, -e[j - 1] / e[1], work, A, n * (j - 1) + 2);
              A_data = A->data;
            }
          }
        }
        for (j = qp1; j < 4; j++) {
          Vf[j - 1] = e[j - 1];
        }
      }
    }
    if (A->size[0] + 1 >= 3) {
      m = 3;
    } else {
      m = A->size[0] + 1;
    }
    if (nct < 3) {
      s_data[nct] = A_data[nct + A->size[0] * nct];
    }
    if (A->size[0] < m) {
      s_data[m - 1] = 0.0;
    }
    if (nrt + 1 < m) {
      e[nrt] = A_data[nrt + A->size[0] * (m - 1)];
    }
    e[m - 1] = 0.0;
    if (nct + 1 <= S_size) {
      for (j = nctp1; j <= S_size; j++) {
        for (k = 0; k < n; k++) {
          U_data[k + U->size[0] * (j - 1)] = 0.0;
        }
        U_data[(j + U->size[0] * (j - 1)) - 1] = 1.0;
      }
    }
    for (q = nct; q >= 1; q--) {
      qp1 = q + 1;
      nmq = n - q;
      qq = (q + n * (q - 1)) - 1;
      if (s_data[q - 1] != 0.0) {
        for (j = qp1; j <= S_size; j++) {
          ns = q + n * (j - 1);
          nrm = 0.0;
          if (nmq >= 0) {
            for (k = 0; k <= nmq; k++) {
              nrm += U_data[qq + k] * U_data[(ns + k) - 1];
            }
          }
          nrm = -(nrm / U_data[qq]);
          xaxpy(nmq + 1, nrm, qq + 1, U, ns);
          U_data = U->data;
        }
        ns = (((nmq + 1) / 2) << 1) + q;
        iy = ns - 2;
        for (j = q; j <= iy; j += 2) {
          r = _mm_loadu_pd(&U_data[(j + U->size[0] * (q - 1)) - 1]);
          _mm_storeu_pd(&U_data[(j + U->size[0] * (q - 1)) - 1],
                        _mm_mul_pd(r, _mm_set1_pd(-1.0)));
        }
        for (j = ns; j <= n; j++) {
          U_data[(j + U->size[0] * (q - 1)) - 1] =
              -U_data[(j + U->size[0] * (q - 1)) - 1];
        }
        U_data[qq]++;
        for (j = 0; j <= q - 2; j++) {
          U_data[j + U->size[0] * (q - 1)] = 0.0;
        }
      } else {
        for (j = 0; j < n; j++) {
          U_data[j + U->size[0] * (q - 1)] = 0.0;
        }
        U_data[qq] = 1.0;
      }
    }
    for (j = 2; j >= 0; j--) {
      if ((j + 1 <= nrt) && (e[0] != 0.0)) {
        b_xaxpy(-(xdotc(Vf, Vf, 5) / Vf[1]), Vf, 5);
        b_xaxpy(-(xdotc(Vf, Vf, 8) / Vf[1]), Vf, 8);
      }
      Vf[3 * j] = 0.0;
      Vf[3 * j + 1] = 0.0;
      Vf[3 * j + 2] = 0.0;
      Vf[j + 3 * j] = 1.0;
    }
    nctp1 = (unsigned char)m;
    for (k = 0; k < nctp1; k++) {
      nrm = s_data[k];
      if (nrm != 0.0) {
        rt = fabs(nrm);
        nrm /= rt;
        s_data[k] = rt;
        if (k + 1 < m) {
          e[k] /= nrm;
        }
        if (k + 1 <= n) {
          ns = n * k;
          iy = ns + n;
          qq = ((((iy - ns) / 2) << 1) + ns) + 1;
          ix = qq - 2;
          for (j = ns + 1; j <= ix; j += 2) {
            r = _mm_loadu_pd(&U_data[j - 1]);
            _mm_storeu_pd(&U_data[j - 1], _mm_mul_pd(_mm_set1_pd(nrm), r));
          }
          for (j = qq; j <= iy; j++) {
            U_data[j - 1] *= nrm;
          }
        }
      }
      if (k + 1 < m) {
        nrm = e[k];
        if (nrm != 0.0) {
          rt = fabs(nrm);
          nrm = rt / nrm;
          e[k] = rt;
          s_data[k + 1] *= nrm;
          ns = 3 * (k + 1);
          iy = ns + 3;
          qq = ns + 3;
          ix = ns + 1;
          for (j = ns + 1; j <= ix; j += 2) {
            r = _mm_loadu_pd(&Vf[j - 1]);
            _mm_storeu_pd(&Vf[j - 1], _mm_mul_pd(_mm_set1_pd(nrm), r));
          }
          for (j = qq; j <= iy; j++) {
            Vf[j - 1] *= nrm;
          }
        }
      }
    }
    nmq = m;
    iter = 0;
    snorm = 0.0;
    for (j = 0; j < nctp1; j++) {
      snorm = fmax(snorm, fmax(fabs(s_data[j]), fabs(e[j])));
    }
    while ((m > 0) && (iter < 75)) {
      boolean_T exitg1;
      nctp1 = m - 1;
      exitg1 = false;
      while (!(exitg1 || (nctp1 == 0))) {
        nrm = fabs(e[nctp1 - 1]);
        if ((nrm <= 2.2204460492503131E-16 *
                        (fabs(s_data[nctp1 - 1]) + fabs(s_data[nctp1]))) ||
            (nrm <= 1.0020841800044864E-292) ||
            ((iter > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
          e[nctp1 - 1] = 0.0;
          exitg1 = true;
        } else {
          nctp1--;
        }
      }
      if (nctp1 == m - 1) {
        ns = 4;
      } else {
        iy = m;
        ns = m;
        exitg1 = false;
        while ((!exitg1) && (ns >= nctp1)) {
          iy = ns;
          if (ns == nctp1) {
            exitg1 = true;
          } else {
            nrm = 0.0;
            if (ns < m) {
              nrm = fabs(e[ns - 1]);
            }
            if (ns > nctp1 + 1) {
              nrm += fabs(e[ns - 2]);
            }
            rt = fabs(s_data[ns - 1]);
            if ((rt <= 2.2204460492503131E-16 * nrm) ||
                (rt <= 1.0020841800044864E-292)) {
              s_data[ns - 1] = 0.0;
              exitg1 = true;
            } else {
              ns--;
            }
          }
        }
        if (iy == nctp1) {
          ns = 3;
        } else if (iy == m) {
          ns = 1;
        } else {
          ns = 2;
          nctp1 = iy;
        }
      }
      switch (ns) {
      case 1:
        f = e[m - 2];
        e[m - 2] = 0.0;
        ns = m - 1;
        for (j = ns; j >= nctp1 + 1; j--) {
          rt = xrotg(&s_data[j - 1], &f, &nrm);
          if (j > nctp1 + 1) {
            f = -nrm * e[0];
            e[0] *= rt;
          }
          xrot(Vf, 3 * (j - 1) + 1, 3 * (m - 1) + 1, rt, nrm);
        }
        break;
      case 2:
        f = e[nctp1 - 1];
        e[nctp1 - 1] = 0.0;
        for (j = nctp1 + 1; j <= m; j++) {
          sqds = xrotg(&s_data[j - 1], &f, &b);
          nrm = e[j - 1];
          f = -b * nrm;
          e[j - 1] = nrm * sqds;
          if (n >= 1) {
            ns = n * (j - 1);
            qq = n * (nctp1 - 1);
            for (k = 0; k < n; k++) {
              ix = qq + k;
              nrm = U_data[ix];
              iy = ns + k;
              rt = U_data[iy];
              U_data[ix] = sqds * nrm - b * rt;
              U_data[iy] = sqds * rt + b * nrm;
            }
          }
        }
        break;
      case 3: {
        double scale;
        nrm = s_data[m - 1];
        rt = s_data[m - 2];
        b = e[m - 2];
        scale = fmax(
            fmax(fmax(fmax(fabs(nrm), fabs(rt)), fabs(b)), fabs(s_data[nctp1])),
            fabs(e[nctp1]));
        sm = nrm / scale;
        nrm = rt / scale;
        rt = b / scale;
        sqds = s_data[nctp1] / scale;
        b = ((nrm + sm) * (nrm - sm) + rt * rt) / 2.0;
        nrm = sm * rt;
        nrm *= nrm;
        if ((b != 0.0) || (nrm != 0.0)) {
          rt = sqrt(b * b + nrm);
          if (b < 0.0) {
            rt = -rt;
          }
          rt = nrm / (b + rt);
        } else {
          rt = 0.0;
        }
        f = (sqds + sm) * (sqds - sm) + rt;
        sqds *= e[nctp1] / scale;
        for (j = nctp1 + 1; j < m; j++) {
          b = xrotg(&f, &sqds, &sm);
          if (j > nctp1 + 1) {
            e[0] = f;
          }
          nrm = e[j - 1];
          rt = s_data[j - 1];
          e[j - 1] = b * nrm - sm * rt;
          sqds = sm * s_data[j];
          s_data[j] *= b;
          xrot(Vf, 3 * (j - 1) + 1, 3 * j + 1, b, sm);
          s_data[j - 1] = b * rt + sm * nrm;
          b = xrotg(&s_data[j - 1], &sqds, &sm);
          nrm = e[j - 1];
          f = b * nrm + sm * s_data[j];
          s_data[j] = -sm * nrm + b * s_data[j];
          sqds = sm * e[j];
          e[j] *= b;
          if (j < n) {
            ns = n * (j - 1);
            iy = n * j;
            for (k = 0; k < n; k++) {
              qq = iy + k;
              nrm = U_data[qq];
              ix = ns + k;
              rt = U_data[ix];
              U_data[qq] = b * nrm - sm * rt;
              U_data[ix] = b * rt + sm * nrm;
            }
          }
        }
        e[m - 2] = f;
        iter++;
      } break;
      default:
        if (s_data[nctp1] < 0.0) {
          s_data[nctp1] = -s_data[nctp1];
          ns = 3 * nctp1;
          iy = ns + 3;
          qq = ns + 3;
          ix = ns + 1;
          for (j = ns + 1; j <= ix; j += 2) {
            r = _mm_loadu_pd(&Vf[j - 1]);
            _mm_storeu_pd(&Vf[j - 1], _mm_mul_pd(r, _mm_set1_pd(-1.0)));
          }
          for (j = qq; j <= iy; j++) {
            Vf[j - 1] = -Vf[j - 1];
          }
        }
        qp1 = nctp1 + 1;
        while ((nctp1 + 1 < nmq) && (s_data[nctp1] < s_data[qp1])) {
          rt = s_data[nctp1];
          s_data[nctp1] = s_data[qp1];
          s_data[qp1] = rt;
          xswap(Vf, 3 * nctp1 + 1, 3 * (nctp1 + 1) + 1);
          if (nctp1 + 1 < n) {
            ix = n * nctp1;
            ns = n * (nctp1 + 1);
            for (j = 0; j < n; j++) {
              iy = ix + j;
              nrm = U_data[iy];
              qq = ns + j;
              U_data[iy] = U_data[qq];
              U_data[qq] = nrm;
            }
          }
          nctp1 = qp1;
          qp1++;
        }
        iter = 0;
        m--;
        break;
      }
    }
  }
  emxFree_real_T(&work);
  if (S_size - 1 >= 0) {
    memcpy(&S_data[0], &s_data[0], (unsigned int)S_size * sizeof(double));
  }
  if (doscale) {
    xzlascl(cscale, anrm, S_size, S_data);
  }
  V_size[0] = 3;
  V_size[1] = S_size;
  for (j = 0; j < S_size; j++) {
    V_data[3 * j] = Vf[3 * j];
    ns = 3 * j + 1;
    V_data[ns] = Vf[ns];
    ns = 3 * j + 2;
    V_data[ns] = Vf[ns];
  }
  return S_size;
}

/* End of code generation (xzsvdc.c) */
