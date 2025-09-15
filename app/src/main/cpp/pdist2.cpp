//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// pdist2.cpp
//
// Code generation for function 'pdist2'
//

// Include files
#include "pdist2.h"
#include "rt_nonfinite.h"
#include "scanfornan.h"
#include "topkrows.h"
#include <cmath>

// Function Definitions
namespace coder {
namespace statslib {
namespace coder {
double pdist2(const double Xin[7232], const double Yin[32], int &b_I)
{
  double X[7232];
  double d[226];
  double D;
  bool logIndX[226];
  for (int i{0}; i < 226; i++) {
    for (int jj{0}; jj < 32; jj++) {
      X[jj + (i << 5)] = Xin[i + 226 * jj];
    }
    d[i] = rtNaN;
  }
  pdistutils::scanfornan(X, logIndX);
  if (pdistutils::b_scanfornan(Yin)) {
    for (int i{0}; i < 226; i++) {
      if (logIndX[i]) {
        D = 0.0;
        for (int jj{0}; jj < 32; jj++) {
          double tempSum_tmp;
          tempSum_tmp = X[jj + (i << 5)] - Yin[jj];
          D += tempSum_tmp * tempSum_tmp;
        }
        d[i] = std::sqrt(D);
      }
    }
  }
  return TopKRows(d, b_I);
}

} // namespace coder
} // namespace statslib
} // namespace coder

// End of code generation (pdist2.cpp)
