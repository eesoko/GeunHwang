//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// scanfornan.cpp
//
// Code generation for function 'scanfornan'
//

// Include files
#include "scanfornan.h"
#include "rt_nonfinite.h"
#include "omp.h"
#include <cmath>

// Function Definitions
namespace coder {
namespace statslib {
namespace coder {
namespace pdistutils {
bool b_scanfornan(const double X[32])
{
  int jj;
  bool exitg1;
  bool nanobs;
  nanobs = false;
  jj = 0;
  exitg1 = false;
  while ((!exitg1) && (jj < 32)) {
    if (std::isnan(X[jj])) {
      nanobs = true;
      exitg1 = true;
    } else {
      jj++;
    }
  }
  return !nanobs;
}

void scanfornan(const double X[7232], bool nanobs[226])
{
  int jj;
  bool exitg1;
  bool nanflag;
  for (int i{0}; i < 226; i++) {
    nanobs[i] = true;
  }
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        nanflag, jj, exitg1)

  for (int qq = 0; qq < 226; qq++) {
    nanflag = false;
    jj = 0;
    exitg1 = false;
    while ((!exitg1) && (jj < 32)) {
      if (std::isnan(X[jj + (qq << 5)])) {
        nanflag = true;
        exitg1 = true;
      } else {
        jj++;
      }
    }
    if (nanflag) {
      nanobs[qq] = false;
    }
  }
}

} // namespace pdistutils
} // namespace coder
} // namespace statslib
} // namespace coder

// End of code generation (scanfornan.cpp)
