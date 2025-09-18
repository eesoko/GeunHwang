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
int xzsvdc(double A_data[], const int A_size[2], double U_data[], int U_size[2],
           double S_data[], double V_data[], int V_size[2])
{
  double work_data[150];
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
  int S_size;
  int j;
  int n;
  int ns;
  int q;
  boolean_T doscale;
  n = A_size[0];
  ns = A_size[0] + 1;
  if (ns > 3) {
    ns = 3;
  }
  S_size = A_size[0];
  if (S_size > 3) {
    S_size = 3;
  }
  memset(&s_data[0], 0, (unsigned int)ns * sizeof(double));
  e[0] = 0.0;
  e[1] = 0.0;
  e[2] = 0.0;
  if (n - 1 >= 0) {
    memset(&work_data[0], 0, (unsigned int)n * sizeof(double));
  }
  U_size[0] = A_size[0];
  U_size[1] = S_size;
  ns = A_size[0] * S_size;
  if (ns - 1 >= 0) {
    memset(&U_data[0], 0, (unsigned int)ns * sizeof(double));
  }
  memset(&Vf[0], 0, 9U * sizeof(double));
  doscale = false;
  cscale = 0.0;
  anrm = xzlangeM(A_data, A_size);
  if (A_size[0] == 0) {
    Vf[0] = 1.0;
    Vf[4] = 1.0;
    Vf[8] = 1.0;
  } else {
    __m128d r;
    double rt;
    double snorm;
    int b_vectorUB;
    int iter;
    int m;
    int nct;
    int nctp1;
    int nmq;
    int nrt;
    int qp1;
    int qq;
    int qq_tmp;
    int qs;
    int vectorUB;
    cscale = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      doscale = true;
      cscale = 6.7178761075670888E-139;
      b_xzlascl(anrm, cscale, A_size[0], A_data, A_size[0]);
    } else if (anrm > 1.4885657073574029E+138) {
      doscale = true;
      cscale = 1.4885657073574029E+138;
      b_xzlascl(anrm, cscale, A_size[0], A_data, A_size[0]);
    }
    nrt = (A_size[0] >= 1);
    if (A_size[0] >= 1) {
      nct = A_size[0] - 1;
    } else {
      nct = 0;
    }
    if (nct > 3) {
      nct = 3;
    }
    nctp1 = nct + 1;
    if (nct >= nrt) {
      m = nct;
    } else {
      m = 1;
    }
    for (q = 0; q < m; q++) {
      boolean_T apply_transform;
      qp1 = q + 2;
      iter = n * q;
      qq_tmp = q + iter;
      qq = qq_tmp + 1;
      nmq = n - q;
      apply_transform = false;
      if (q + 1 <= nct) {
        nrm = xnrm2(nmq, A_data, qq_tmp + 1);
        if (nrm > 0.0) {
          apply_transform = true;
          if (A_data[qq_tmp] < 0.0) {
            nrm = -nrm;
          }
          s_data[q] = nrm;
          if (fabs(nrm) >= 1.0020841800044864E-292) {
            nrm = 1.0 / nrm;
            ns = (qq_tmp + nmq) + 1;
            qs = (((((ns - qq_tmp) - 1) / 2) << 1) + qq_tmp) + 1;
            vectorUB = qs - 2;
            for (j = qq; j <= vectorUB; j += 2) {
              r = _mm_loadu_pd(&A_data[j - 1]);
              _mm_storeu_pd(&A_data[j - 1], _mm_mul_pd(_mm_set1_pd(nrm), r));
            }
            for (j = qs; j < ns; j++) {
              A_data[j - 1] *= nrm;
            }
          } else {
            ns = (qq_tmp + nmq) + 1;
            qs = (((((ns - qq_tmp) - 1) / 2) << 1) + qq_tmp) + 1;
            b_vectorUB = qs - 2;
            for (j = qq; j <= b_vectorUB; j += 2) {
              r = _mm_loadu_pd(&A_data[j - 1]);
              _mm_storeu_pd(&A_data[j - 1],
                            _mm_div_pd(r, _mm_set1_pd(s_data[q])));
            }
            for (j = qs; j < ns; j++) {
              A_data[j - 1] /= s_data[q];
            }
          }
          A_data[qq_tmp]++;
          s_data[q] = -s_data[q];
        } else {
          s_data[q] = 0.0;
        }
      }
      for (j = qp1; j < 4; j++) {
        ns = q + n * (j - 1);
        if (apply_transform) {
          xaxpy(nmq,
                -(xdotc(nmq, A_data, qq_tmp + 1, A_data, ns + 1) /
                  A_data[qq_tmp]),
                qq_tmp + 1, A_data, ns + 1);
        }
        e[j - 1] = A_data[ns];
      }
      if (q + 1 <= nct) {
        for (j = q + 1; j <= n; j++) {
          ns = (j + iter) - 1;
          U_data[ns] = A_data[ns];
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
            qs = ns - 2;
            for (j = qp1; j <= qs; j += 2) {
              r = _mm_loadu_pd(&e[j - 1]);
              _mm_storeu_pd(&e[j - 1], _mm_mul_pd(_mm_set1_pd(nrm), r));
            }
            for (j = ns; j < 4; j++) {
              e[j - 1] *= nrm;
            }
          } else {
            ns = ((((2 - q) / 2) << 1) + q) + 2;
            qs = ns - 2;
            for (j = qp1; j <= qs; j += 2) {
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
            if (qp1 <= n) {
              memset(&work_data[qp1 + -1], 0,
                     (unsigned int)((n - qp1) + 1) * sizeof(double));
            }
            for (j = qp1; j < 4; j++) {
              c_xaxpy(nmq - 1, e[j - 1], A_data, n * (j - 1) + 2, work_data);
            }
            for (j = qp1; j < 4; j++) {
              d_xaxpy(nmq - 1, -e[j - 1] / e[1], work_data, A_data,
                      n * (j - 1) + 2);
            }
          }
        }
        for (j = qp1; j < 4; j++) {
          Vf[j - 1] = e[j - 1];
        }
      }
    }
    if (A_size[0] + 1 >= 3) {
      m = 2;
    } else {
      m = A_size[0];
    }
    if (nct < 3) {
      s_data[nct] = A_data[nct + A_size[0] * nct];
    }
    if (A_size[0] < m + 1) {
      s_data[m] = 0.0;
    }
    if (nrt < m) {
      e[nrt] = A_data[nrt + A_size[0] * m];
    }
    e[m] = 0.0;
    if (nct + 1 <= S_size) {
      for (j = nctp1; j <= S_size; j++) {
        for (q = 0; q < n; q++) {
          U_data[q + U_size[0] * (j - 1)] = 0.0;
        }
        U_data[(j + U_size[0] * (j - 1)) - 1] = 1.0;
      }
    }
    for (q = nct; q >= 1; q--) {
      qp1 = q + 1;
      nmq = n - q;
      vectorUB = n * (q - 1);
      qq = (q + vectorUB) - 1;
      if (s_data[q - 1] != 0.0) {
        for (j = qp1; j <= S_size; j++) {
          ns = q + n * (j - 1);
          xaxpy(nmq + 1,
                -(xdotc(nmq + 1, U_data, qq + 1, U_data, ns) / U_data[qq]),
                qq + 1, U_data, ns);
        }
        b_vectorUB = (((nmq + 1) / 2) << 1) + q;
        ns = b_vectorUB - 2;
        for (j = q; j <= ns; j += 2) {
          qs = (j + vectorUB) - 1;
          r = _mm_loadu_pd(&U_data[qs]);
          _mm_storeu_pd(&U_data[qs], _mm_mul_pd(r, _mm_set1_pd(-1.0)));
        }
        for (j = b_vectorUB; j <= n; j++) {
          ns = (j + vectorUB) - 1;
          U_data[ns] = -U_data[ns];
        }
        U_data[qq]++;
        for (j = 0; j <= q - 2; j++) {
          U_data[j + vectorUB] = 0.0;
        }
      } else {
        if (n - 1 >= 0) {
          memset(&U_data[vectorUB], 0, (unsigned int)n * sizeof(double));
        }
        U_data[qq] = 1.0;
      }
    }
    for (j = 2; j >= 0; j--) {
      if ((j + 1 <= nrt) && (e[0] != 0.0)) {
        b_xaxpy(-(b_xdotc(Vf, Vf, 5) / Vf[1]), Vf, 5);
        b_xaxpy(-(b_xdotc(Vf, Vf, 8) / Vf[1]), Vf, 8);
      }
      Vf[3 * j] = 0.0;
      Vf[3 * j + 1] = 0.0;
      Vf[3 * j + 2] = 0.0;
      Vf[j + 3 * j] = 1.0;
    }
    qq = m;
    iter = 0;
    snorm = 0.0;
    for (q = 0; q <= m; q++) {
      nrm = s_data[q];
      if (nrm != 0.0) {
        rt = fabs(nrm);
        nrm /= rt;
        s_data[q] = rt;
        if (q < m) {
          e[q] /= nrm;
        }
        if (q + 1 <= n) {
          ns = n * q;
          qs = ns + n;
          b_vectorUB = ((((qs - ns) / 2) << 1) + ns) + 1;
          vectorUB = b_vectorUB - 2;
          for (j = ns + 1; j <= vectorUB; j += 2) {
            r = _mm_loadu_pd(&U_data[j - 1]);
            _mm_storeu_pd(&U_data[j - 1], _mm_mul_pd(_mm_set1_pd(nrm), r));
          }
          for (j = b_vectorUB; j <= qs; j++) {
            U_data[j - 1] *= nrm;
          }
        }
      }
      if (q < m) {
        nrm = e[q];
        if (nrm != 0.0) {
          rt = fabs(nrm);
          nrm = rt / nrm;
          e[q] = rt;
          s_data[q + 1] *= nrm;
          ns = 3 * (q + 1);
          qs = ns + 3;
          b_vectorUB = ns + 3;
          vectorUB = ns + 1;
          for (j = ns + 1; j <= vectorUB; j += 2) {
            r = _mm_loadu_pd(&Vf[j - 1]);
            _mm_storeu_pd(&Vf[j - 1], _mm_mul_pd(_mm_set1_pd(nrm), r));
          }
          for (j = b_vectorUB; j <= qs; j++) {
            Vf[j - 1] *= nrm;
          }
        }
      }
      snorm = fmax(snorm, fmax(fabs(s_data[q]), fabs(e[q])));
    }
    while ((m + 1 > 0) && (iter < 75)) {
      boolean_T exitg1;
      qq_tmp = m;
      exitg1 = false;
      while (!(exitg1 || (qq_tmp == 0))) {
        nrm = fabs(e[qq_tmp - 1]);
        if ((nrm <= 2.2204460492503131E-16 *
                        (fabs(s_data[qq_tmp - 1]) + fabs(s_data[qq_tmp]))) ||
            (nrm <= 1.0020841800044864E-292) ||
            ((iter > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
          e[qq_tmp - 1] = 0.0;
          exitg1 = true;
        } else {
          qq_tmp--;
        }
      }
      if (qq_tmp == m) {
        ns = 4;
      } else {
        qs = m + 1;
        ns = m + 1;
        exitg1 = false;
        while ((!exitg1) && (ns >= qq_tmp)) {
          qs = ns;
          if (ns == qq_tmp) {
            exitg1 = true;
          } else {
            nrm = 0.0;
            if (ns < m + 1) {
              nrm = fabs(e[ns - 1]);
            }
            if (ns > qq_tmp + 1) {
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
        if (qs == qq_tmp) {
          ns = 3;
        } else if (qs == m + 1) {
          ns = 1;
        } else {
          ns = 2;
          qq_tmp = qs;
        }
      }
      switch (ns) {
      case 1:
        f = e[m - 1];
        e[m - 1] = 0.0;
        for (j = m; j >= qq_tmp + 1; j--) {
          rt = xrotg(&s_data[j - 1], &f, &nrm);
          if (j > qq_tmp + 1) {
            f = -nrm * e[0];
            e[0] *= rt;
          }
          xrot(Vf, 3 * (j - 1) + 1, 3 * m + 1, rt, nrm);
        }
        break;
      case 2:
        f = e[qq_tmp - 1];
        e[qq_tmp - 1] = 0.0;
        for (j = qq_tmp + 1; j <= m + 1; j++) {
          rt = xrotg(&s_data[j - 1], &f, &nrm);
          b = e[j - 1];
          f = -nrm * b;
          e[j - 1] = b * rt;
          b_xrot(n, U_data, n * (j - 1) + 1, n * (qq_tmp - 1) + 1, rt, nrm);
        }
        break;
      case 3: {
        double scale;
        nrm = s_data[m - 1];
        rt = e[m - 1];
        scale = fmax(fmax(fmax(fmax(fabs(s_data[m]), fabs(nrm)), fabs(rt)),
                          fabs(s_data[qq_tmp])),
                     fabs(e[qq_tmp]));
        sm = s_data[m] / scale;
        nrm /= scale;
        rt /= scale;
        sqds = s_data[qq_tmp] / scale;
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
        nrm = sqds * (e[qq_tmp] / scale);
        for (j = qq_tmp + 1; j <= m; j++) {
          b = xrotg(&f, &nrm, &sm);
          if (j > qq_tmp + 1) {
            e[0] = f;
          }
          nrm = e[j - 1];
          rt = s_data[j - 1];
          e[j - 1] = b * nrm - sm * rt;
          sqds = sm * s_data[j];
          s_data[j] *= b;
          xrot(Vf, 3 * (j - 1) + 1, 3 * j + 1, b, sm);
          s_data[j - 1] = b * rt + sm * nrm;
          rt = xrotg(&s_data[j - 1], &sqds, &b);
          nrm = e[j - 1];
          f = rt * nrm + b * s_data[j];
          s_data[j] = -b * nrm + rt * s_data[j];
          nrm = b * e[j];
          e[j] *= rt;
          if (j < n) {
            b_xrot(n, U_data, n * (j - 1) + 1, n * j + 1, rt, b);
          }
        }
        e[m - 1] = f;
        iter++;
      } break;
      default:
        if (s_data[qq_tmp] < 0.0) {
          s_data[qq_tmp] = -s_data[qq_tmp];
          ns = 3 * qq_tmp;
          qs = ns + 3;
          b_vectorUB = ns + 3;
          vectorUB = ns + 1;
          for (j = ns + 1; j <= vectorUB; j += 2) {
            r = _mm_loadu_pd(&Vf[j - 1]);
            _mm_storeu_pd(&Vf[j - 1], _mm_mul_pd(r, _mm_set1_pd(-1.0)));
          }
          for (j = b_vectorUB; j <= qs; j++) {
            Vf[j - 1] = -Vf[j - 1];
          }
        }
        qp1 = qq_tmp + 1;
        while ((qq_tmp + 1 < qq + 1) && (s_data[qq_tmp] < s_data[qp1])) {
          rt = s_data[qq_tmp];
          s_data[qq_tmp] = s_data[qp1];
          s_data[qp1] = rt;
          xswap(Vf, 3 * qq_tmp + 1, 3 * (qq_tmp + 1) + 1);
          if (qq_tmp + 1 < n) {
            b_xswap(n, U_data, n * qq_tmp + 1, n * (qq_tmp + 1) + 1);
          }
          qq_tmp = qp1;
          qp1++;
        }
        iter = 0;
        m--;
        break;
      }
    }
  }
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
