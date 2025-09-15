//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// xzsvdc.cpp
//
// Code generation for function 'xzsvdc'
//

// Include files
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
#include "coder_array.h"
#include <algorithm>
#include <cmath>
#include <cstring>

// Function Definitions
namespace coder {
namespace internal {
namespace reflapack {
int xzsvdc(array<double, 2U> &A, array<double, 2U> &U, double S_data[],
           double V_data[], int V_size[2])
{
  array<double, 1U> work;
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
  int n;
  int ns;
  bool doscale;
  n = A.size(0);
  ns = A.size(0) + 1;
  if (ns > 3) {
    ns = 3;
  }
  S_size = A.size(0);
  if (S_size > 3) {
    S_size = 3;
  }
  if (ns - 1 >= 0) {
    std::memset(&s_data[0], 0, static_cast<unsigned int>(ns) * sizeof(double));
  }
  e[0] = 0.0;
  e[1] = 0.0;
  e[2] = 0.0;
  work.set_size(A.size(0));
  if (n - 1 >= 0) {
    std::memset(&work[0], 0, static_cast<unsigned int>(n) * sizeof(double));
  }
  U.set_size(A.size(0), S_size);
  ns = A.size(0) * S_size;
  for (int j{0}; j < ns; j++) {
    U[j] = 0.0;
  }
  std::memset(&Vf[0], 0, 9U * sizeof(double));
  doscale = false;
  cscale = 0.0;
  anrm = xzlangeM(A);
  if (A.size(0) == 0) {
    Vf[0] = 1.0;
    Vf[4] = 1.0;
    Vf[8] = 1.0;
  } else {
    double rt;
    double snorm;
    int iy;
    int m;
    int nct;
    int nctp1;
    int nmq;
    int nrt;
    int qp1;
    int qq;
    cscale = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      doscale = true;
      cscale = 6.7178761075670888E-139;
      xzlascl(anrm, cscale, A.size(0), A, A.size(0));
    } else if (anrm > 1.4885657073574029E+138) {
      doscale = true;
      cscale = 1.4885657073574029E+138;
      xzlascl(anrm, cscale, A.size(0), A, A.size(0));
    }
    nrt = (A.size(0) >= 1);
    if (A.size(0) >= 1) {
      nct = A.size(0) - 1;
    } else {
      nct = 0;
    }
    if (nct > 3) {
      nct = 3;
    }
    nctp1 = nct + 1;
    if (nct >= nrt) {
      iy = nct;
    } else {
      iy = 1;
    }
    for (int q{0}; q < iy; q++) {
      bool apply_transform;
      qp1 = q + 2;
      qq = (q + n * q) + 1;
      nmq = (n - q) - 1;
      apply_transform = false;
      if (q + 1 <= nct) {
        nrm = blas::xnrm2(nmq + 1, A, qq);
        if (nrm > 0.0) {
          apply_transform = true;
          if (A[qq - 1] < 0.0) {
            nrm = -nrm;
          }
          s_data[q] = nrm;
          if (std::abs(nrm) >= 1.0020841800044864E-292) {
            nrm = 1.0 / nrm;
            ns = qq + nmq;
            for (int j{qq}; j <= ns; j++) {
              A[j - 1] = nrm * A[j - 1];
            }
          } else {
            ns = qq + nmq;
            for (int j{qq}; j <= ns; j++) {
              A[j - 1] = A[j - 1] / s_data[q];
            }
          }
          A[qq - 1] = A[qq - 1] + 1.0;
          s_data[q] = -s_data[q];
        } else {
          s_data[q] = 0.0;
        }
      }
      for (int j{qp1}; j < 4; j++) {
        ns = q + n * (j - 1);
        if (apply_transform) {
          nrm = 0.0;
          if (nmq >= 0) {
            for (int k{0}; k <= nmq; k++) {
              nrm += A[(qq + k) - 1] * A[ns + k];
            }
          }
          nrm = -(nrm / A[q + A.size(0) * q]);
          blas::xaxpy(nmq + 1, nrm, qq, A, ns + 1);
        }
        e[j - 1] = A[ns];
      }
      if (q + 1 <= nct) {
        for (int j{q + 1}; j <= n; j++) {
          U[(j + U.size(0) * q) - 1] = A[(j + A.size(0) * q) - 1];
        }
      }
      if (q + 1 <= nrt) {
        nrm = blas::xnrm2(e);
        if (nrm == 0.0) {
          e[0] = 0.0;
        } else {
          if (e[1] < 0.0) {
            e[0] = -nrm;
          } else {
            e[0] = nrm;
          }
          nrm = e[0];
          if (std::abs(e[0]) >= 1.0020841800044864E-292) {
            nrm = 1.0 / e[0];
            for (int j{qp1}; j < 4; j++) {
              e[j - 1] *= nrm;
            }
          } else {
            for (int j{qp1}; j < 4; j++) {
              e[j - 1] /= nrm;
            }
          }
          e[1]++;
          e[0] = -e[0];
          if (n >= 2) {
            if (qp1 <= n) {
              std::memset(&work[qp1 + -1], 0,
                          static_cast<unsigned int>((n - qp1) + 1) *
                              sizeof(double));
            }
            for (int j{qp1}; j < 4; j++) {
              blas::xaxpy(nmq, e[j - 1], A, n * (j - 1) + 2, work);
            }
            for (int j{qp1}; j < 4; j++) {
              blas::xaxpy(nmq, -e[j - 1] / e[1], work, A, n * (j - 1) + 2);
            }
          }
        }
        for (int j{qp1}; j < 4; j++) {
          Vf[j - 1] = e[j - 1];
        }
      }
    }
    if (A.size(0) + 1 >= 3) {
      m = 3;
    } else {
      m = A.size(0) + 1;
    }
    if (nct < 3) {
      s_data[nct] = A[nct + A.size(0) * nct];
    }
    if (A.size(0) < m) {
      s_data[m - 1] = 0.0;
    }
    if (nrt + 1 < m) {
      e[nrt] = A[nrt + A.size(0) * (m - 1)];
    }
    e[m - 1] = 0.0;
    if (nct + 1 <= S_size) {
      for (int j{nctp1}; j <= S_size; j++) {
        for (int k{0}; k < n; k++) {
          U[k + U.size(0) * (j - 1)] = 0.0;
        }
        U[(j + U.size(0) * (j - 1)) - 1] = 1.0;
      }
    }
    for (int q{nct}; q >= 1; q--) {
      qp1 = q + 1;
      nmq = n - q;
      qq = (q + n * (q - 1)) - 1;
      if (s_data[q - 1] != 0.0) {
        for (int j{qp1}; j <= S_size; j++) {
          ns = q + n * (j - 1);
          nrm = 0.0;
          if (nmq >= 0) {
            for (int k{0}; k <= nmq; k++) {
              nrm += U[qq + k] * U[(ns + k) - 1];
            }
          }
          nrm = -(nrm / U[qq]);
          blas::xaxpy(nmq + 1, nrm, qq + 1, U, ns);
        }
        for (int j{q}; j <= n; j++) {
          ns = (j + U.size(0) * (q - 1)) - 1;
          U[ns] = -U[ns];
        }
        U[qq] = U[qq] + 1.0;
        for (int j{0}; j <= q - 2; j++) {
          U[j + U.size(0) * (q - 1)] = 0.0;
        }
      } else {
        for (int j{0}; j < n; j++) {
          U[j + U.size(0) * (q - 1)] = 0.0;
        }
        U[qq] = 1.0;
      }
    }
    for (int j{2}; j >= 0; j--) {
      if ((j + 1 <= nrt) && (e[0] != 0.0)) {
        blas::xaxpy(-(blas::xdotc(Vf, Vf, 5) / Vf[1]), Vf, 5);
        blas::xaxpy(-(blas::xdotc(Vf, Vf, 8) / Vf[1]), Vf, 8);
      }
      Vf[3 * j] = 0.0;
      Vf[3 * j + 1] = 0.0;
      Vf[3 * j + 2] = 0.0;
      Vf[j + 3 * j] = 1.0;
    }
    qq = static_cast<unsigned char>(m);
    nrt = m;
    nct = 0;
    snorm = 0.0;
    for (int k{0}; k < qq; k++) {
      nrm = s_data[k];
      if (nrm != 0.0) {
        rt = std::abs(nrm);
        nrm /= rt;
        s_data[k] = rt;
        if (k + 1 < m) {
          e[k] /= nrm;
        }
        if (k + 1 <= n) {
          ns = n * k;
          iy = ns + n;
          for (int j{ns + 1}; j <= iy; j++) {
            U[j - 1] = nrm * U[j - 1];
          }
        }
      }
      if (k + 1 < m) {
        nrm = e[k];
        if (nrm != 0.0) {
          rt = std::abs(nrm);
          nrm = rt / nrm;
          e[k] = rt;
          s_data[k + 1] *= nrm;
          ns = 3 * (k + 1);
          iy = ns + 3;
          for (int j{ns + 1}; j <= iy; j++) {
            Vf[j - 1] *= nrm;
          }
        }
      }
      snorm = std::fmax(snorm, std::fmax(std::abs(s_data[k]), std::abs(e[k])));
    }
    while ((m > 0) && (nct < 75)) {
      bool exitg1;
      nmq = m - 1;
      exitg1 = false;
      while (!(exitg1 || (nmq == 0))) {
        nrm = std::abs(e[nmq - 1]);
        if ((nrm <= 2.2204460492503131E-16 *
                        (std::abs(s_data[nmq - 1]) + std::abs(s_data[nmq]))) ||
            (nrm <= 1.0020841800044864E-292) ||
            ((nct > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
          e[nmq - 1] = 0.0;
          exitg1 = true;
        } else {
          nmq--;
        }
      }
      if (nmq == m - 1) {
        ns = 4;
      } else {
        iy = m;
        ns = m;
        exitg1 = false;
        while ((!exitg1) && (ns >= nmq)) {
          iy = ns;
          if (ns == nmq) {
            exitg1 = true;
          } else {
            nrm = 0.0;
            if (ns < m) {
              nrm = std::abs(e[ns - 1]);
            }
            if (ns > nmq + 1) {
              nrm += std::abs(e[ns - 2]);
            }
            rt = std::abs(s_data[ns - 1]);
            if ((rt <= 2.2204460492503131E-16 * nrm) ||
                (rt <= 1.0020841800044864E-292)) {
              s_data[ns - 1] = 0.0;
              exitg1 = true;
            } else {
              ns--;
            }
          }
        }
        if (iy == nmq) {
          ns = 3;
        } else if (iy == m) {
          ns = 1;
        } else {
          ns = 2;
          nmq = iy;
        }
      }
      switch (ns) {
      case 1:
        f = e[m - 2];
        e[m - 2] = 0.0;
        ns = m - 1;
        for (int j{ns}; j >= nmq + 1; j--) {
          rt = blas::xrotg(s_data[j - 1], f, nrm);
          if (j > nmq + 1) {
            f = -nrm * e[0];
            e[0] *= rt;
          }
          blas::xrot(Vf, 3 * (j - 1) + 1, 3 * (m - 1) + 1, rt, nrm);
        }
        break;
      case 2:
        f = e[nmq - 1];
        e[nmq - 1] = 0.0;
        for (int j{nmq + 1}; j <= m; j++) {
          sqds = blas::xrotg(s_data[j - 1], f, b);
          nrm = e[j - 1];
          f = -b * nrm;
          e[j - 1] = nrm * sqds;
          if (n >= 1) {
            ns = n * (j - 1);
            nctp1 = n * (nmq - 1);
            for (int k{0}; k < n; k++) {
              iy = nctp1 + k;
              nrm = U[iy];
              qq = ns + k;
              rt = U[qq];
              U[iy] = sqds * nrm - b * rt;
              U[qq] = sqds * rt + b * nrm;
            }
          }
        }
        break;
      case 3: {
        double scale;
        nrm = s_data[m - 1];
        rt = s_data[m - 2];
        b = e[m - 2];
        scale = std::fmax(
            std::fmax(
                std::fmax(std::fmax(std::abs(nrm), std::abs(rt)), std::abs(b)),
                std::abs(s_data[nmq])),
            std::abs(e[nmq]));
        sm = nrm / scale;
        nrm = rt / scale;
        rt = b / scale;
        sqds = s_data[nmq] / scale;
        b = ((nrm + sm) * (nrm - sm) + rt * rt) / 2.0;
        nrm = sm * rt;
        nrm *= nrm;
        if ((b != 0.0) || (nrm != 0.0)) {
          rt = std::sqrt(b * b + nrm);
          if (b < 0.0) {
            rt = -rt;
          }
          rt = nrm / (b + rt);
        } else {
          rt = 0.0;
        }
        f = (sqds + sm) * (sqds - sm) + rt;
        sqds *= e[nmq] / scale;
        for (int j{nmq + 1}; j < m; j++) {
          b = blas::xrotg(f, sqds, sm);
          if (j > nmq + 1) {
            e[0] = f;
          }
          nrm = e[j - 1];
          rt = s_data[j - 1];
          e[j - 1] = b * nrm - sm * rt;
          sqds = sm * s_data[j];
          s_data[j] *= b;
          blas::xrot(Vf, 3 * (j - 1) + 1, 3 * j + 1, b, sm);
          s_data[j - 1] = b * rt + sm * nrm;
          b = blas::xrotg(s_data[j - 1], sqds, sm);
          nrm = e[j - 1];
          f = b * nrm + sm * s_data[j];
          s_data[j] = -sm * nrm + b * s_data[j];
          sqds = sm * e[j];
          e[j] *= b;
          if (j < n) {
            ns = n * (j - 1);
            iy = n * j;
            for (int k{0}; k < n; k++) {
              qq = iy + k;
              nrm = U[qq];
              nctp1 = ns + k;
              rt = U[nctp1];
              U[qq] = b * nrm - sm * rt;
              U[nctp1] = b * rt + sm * nrm;
            }
          }
        }
        e[m - 2] = f;
        nct++;
      } break;
      default:
        if (s_data[nmq] < 0.0) {
          s_data[nmq] = -s_data[nmq];
          ns = 3 * nmq;
          iy = ns + 3;
          for (int j{ns + 1}; j <= iy; j++) {
            Vf[j - 1] = -Vf[j - 1];
          }
        }
        qp1 = nmq + 1;
        while ((nmq + 1 < nrt) && (s_data[nmq] < s_data[qp1])) {
          rt = s_data[nmq];
          s_data[nmq] = s_data[qp1];
          s_data[qp1] = rt;
          blas::xswap(Vf, 3 * nmq + 1, 3 * (nmq + 1) + 1);
          if (nmq + 1 < n) {
            nctp1 = n * nmq;
            ns = n * (nmq + 1);
            for (int j{0}; j < n; j++) {
              iy = nctp1 + j;
              nrm = U[iy];
              qq = ns + j;
              U[iy] = U[qq];
              U[qq] = nrm;
            }
          }
          nmq = qp1;
          qp1++;
        }
        nct = 0;
        m--;
        break;
      }
    }
  }
  if (S_size - 1 >= 0) {
    std::copy(&s_data[0], &s_data[S_size], &S_data[0]);
  }
  if (doscale) {
    xzlascl(cscale, anrm, S_size, S_data);
  }
  V_size[0] = 3;
  V_size[1] = S_size;
  for (int j{0}; j < S_size; j++) {
    V_data[3 * j] = Vf[3 * j];
    ns = 3 * j + 1;
    V_data[ns] = Vf[ns];
    ns = 3 * j + 2;
    V_data[ns] = Vf[ns];
  }
  return S_size;
}

} // namespace reflapack
} // namespace internal
} // namespace coder

// End of code generation (xzsvdc.cpp)
