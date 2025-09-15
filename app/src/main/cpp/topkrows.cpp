//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// topkrows.cpp
//
// Code generation for function 'topkrows'
//

// Include files
#include "topkrows.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions
namespace coder {
double TopKRows(const double X[226], int &b_I)
{
  double B;
  int heap;
  heap = 0;
  for (int ii{0}; ii < 225; ii++) {
    bool newRowIsHigher;
    if (ii + 2 < heap + 1) {
      bool isnanv2;
      newRowIsHigher = true;
      isnanv2 = std::isnan(X[heap]);
      B = X[ii + 1];
      if ((!(B == X[heap])) && ((!std::isnan(B)) || (!isnanv2)) &&
          (!(B <= X[heap])) && (!isnanv2)) {
        newRowIsHigher = false;
      }
    } else {
      bool isnanv2;
      newRowIsHigher = true;
      B = X[ii + 1];
      isnanv2 = std::isnan(B);
      if ((!(X[heap] == B)) && ((!std::isnan(X[heap])) || (!isnanv2)) &&
          (!(X[heap] <= B)) && (!isnanv2)) {
        newRowIsHigher = false;
      }
      newRowIsHigher = !newRowIsHigher;
    }
    if (newRowIsHigher) {
      heap = ii + 1;
    }
  }
  b_I = heap + 1;
  return X[heap];
}

} // namespace coder

// End of code generation (topkrows.cpp)
