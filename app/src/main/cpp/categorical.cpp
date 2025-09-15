//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// categorical.cpp
//
// Code generation for function 'categorical'
//

// Include files
#include "categorical.h"
#include "anonymous_function.h"
#include "find.h"
#include "introsort.h"
#include "predictExercise_internal_types.h"
#include "predictExercise_types.h"
#include "rt_nonfinite.h"
#include "strcmp.h"
#include "coder_array.h"
#include "coder_bounded_array.h"
#include <cstring>

// Function Definitions
namespace coder {
categorical::categorical() = default;

categorical::~categorical() = default;

void categorical::convertCodes(const cell_wrap_3 bnames[6],
                               const cell_wrap_3 anames[6],
                               unsigned char bcodesout[6])
{
  anonymous_function b_this;
  array<signed char, 2U> y;
  cell_wrap_3 c[6];
  cell_wrap_3 uA_data[6];
  cell_wrap_3 uB_data[6];
  double ib_data[6];
  double icA[6];
  double y_data[5];
  double b_d;
  int a__2_data[6];
  int b_ib_data[6];
  int dIdx_data[6];
  int idx[6];
  int b_nz;
  int c_i;
  int d_i;
  int j;
  int loop_ub;
  int minSize;
  int nz;
  unsigned char b2a[7];
  signed char ii_data[6];
  bool d[6];
  bool uAInUB_data[6];
  bool exitg1;
  bool is_less_than;
  for (int i{0}; i < 6; i++) {
    b_this.workspace.c[i] = bnames[i];
  }
  internal::introsort(b_this, idx);
  loop_ub = bnames[idx[0] - 1].f1.size[1];
  c[0].f1.size[0] = 1;
  c[0].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[0].f1.data[i] = bnames[idx[0] - 1].f1.data[i];
  }
  loop_ub = bnames[idx[1] - 1].f1.size[1];
  c[1].f1.size[0] = 1;
  c[1].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[1].f1.data[i] = bnames[idx[1] - 1].f1.data[i];
  }
  loop_ub = bnames[idx[2] - 1].f1.size[1];
  c[2].f1.size[0] = 1;
  c[2].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[2].f1.data[i] = bnames[idx[2] - 1].f1.data[i];
  }
  loop_ub = bnames[idx[3] - 1].f1.size[1];
  c[3].f1.size[0] = 1;
  c[3].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[3].f1.data[i] = bnames[idx[3] - 1].f1.data[i];
  }
  loop_ub = bnames[idx[4] - 1].f1.size[1];
  c[4].f1.size[0] = 1;
  c[4].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[4].f1.data[i] = bnames[idx[4] - 1].f1.data[i];
  }
  loop_ub = bnames[idx[5] - 1].f1.size[1];
  c[5].f1.size[0] = 1;
  c[5].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[5].f1.data[i] = bnames[idx[5] - 1].f1.data[i];
  }
  d[0] = true;
  nz = 1;
  for (int i{0}; i < 5; i++) {
    is_less_than = !internal::b_strcmp(c[i].f1.data, c[i].f1.size,
                                       c[i + 1].f1.data, c[i + 1].f1.size);
    d[i + 1] = is_less_than;
    nz += is_less_than;
  }
  j = eml_find(d, dIdx_data);
  for (int b_i{0}; b_i < j; b_i++) {
    if (b_i + 1 != j) {
      b_d = static_cast<double>(dIdx_data[b_i + 1]) - 1.0;
      loop_ub = dIdx_data[b_i];
      if (b_d < loop_ub) {
        c_i = 0;
      } else {
        minSize = static_cast<int>(b_d) - loop_ub;
        c_i = minSize + 1;
        for (int i{0}; i <= minSize; i++) {
          y_data[i] = static_cast<double>(loop_ub) + static_cast<double>(i);
        }
      }
      loop_ub = idx[static_cast<int>(y_data[0]) - 1];
      for (int i{2}; i <= c_i; i++) {
        minSize = idx[static_cast<int>(y_data[i - 1]) - 1];
        if (loop_ub > minSize) {
          loop_ub = minSize;
        }
      }
      a__2_data[b_i] = loop_ub;
    } else {
      loop_ub = dIdx_data[b_i];
      c_i = 7 - loop_ub;
      loop_ub = 6 - loop_ub;
      for (int i{0}; i <= loop_ub; i++) {
        ib_data[i] =
            static_cast<double>(dIdx_data[b_i]) + static_cast<double>(i);
      }
      loop_ub = idx[static_cast<int>(ib_data[0]) - 1];
      for (int i{2}; i <= c_i; i++) {
        minSize = idx[static_cast<int>(ib_data[i - 1]) - 1];
        if (loop_ub > minSize) {
          loop_ub = minSize;
        }
      }
      a__2_data[b_i] = loop_ub;
    }
  }
  for (int i{0}; i < 6; i++) {
    icA[i] = d[i];
  }
  for (int i{0}; i < 5; i++) {
    icA[i + 1] += icA[i];
  }
  for (int i{0}; i < 6; i++) {
    ib_data[i] = icA[i];
  }
  for (int i{0}; i < 6; i++) {
    icA[idx[i] - 1] = ib_data[i];
  }
  for (int i{0}; i < j; i++) {
    loop_ub = bnames[a__2_data[i] - 1].f1.size[1];
    uA_data[i].f1.size[0] = 1;
    uA_data[i].f1.size[1] = loop_ub;
    for (int b_i{0}; b_i < loop_ub; b_i++) {
      uA_data[i].f1.data[b_i] = bnames[a__2_data[i] - 1].f1.data[b_i];
    }
  }
  for (int i{0}; i < 6; i++) {
    b_this.workspace.c[i] = anames[i];
  }
  internal::introsort(b_this, idx);
  loop_ub = anames[idx[0] - 1].f1.size[1];
  c[0].f1.size[0] = 1;
  c[0].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[0].f1.data[i] = anames[idx[0] - 1].f1.data[i];
  }
  loop_ub = anames[idx[1] - 1].f1.size[1];
  c[1].f1.size[0] = 1;
  c[1].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[1].f1.data[i] = anames[idx[1] - 1].f1.data[i];
  }
  loop_ub = anames[idx[2] - 1].f1.size[1];
  c[2].f1.size[0] = 1;
  c[2].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[2].f1.data[i] = anames[idx[2] - 1].f1.data[i];
  }
  loop_ub = anames[idx[3] - 1].f1.size[1];
  c[3].f1.size[0] = 1;
  c[3].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[3].f1.data[i] = anames[idx[3] - 1].f1.data[i];
  }
  loop_ub = anames[idx[4] - 1].f1.size[1];
  c[4].f1.size[0] = 1;
  c[4].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[4].f1.data[i] = anames[idx[4] - 1].f1.data[i];
  }
  loop_ub = anames[idx[5] - 1].f1.size[1];
  c[5].f1.size[0] = 1;
  c[5].f1.size[1] = loop_ub;
  for (int i{0}; i < loop_ub; i++) {
    c[5].f1.data[i] = anames[idx[5] - 1].f1.data[i];
  }
  d[0] = true;
  b_nz = 1;
  for (int i{0}; i < 5; i++) {
    is_less_than = !internal::b_strcmp(c[i].f1.data, c[i].f1.size,
                                       c[i + 1].f1.data, c[i + 1].f1.size);
    d[i + 1] = is_less_than;
    b_nz += is_less_than;
  }
  j = eml_find(d, dIdx_data);
  for (int b_i{0}; b_i < j; b_i++) {
    if (b_i + 1 != j) {
      b_d = static_cast<double>(dIdx_data[b_i + 1]) - 1.0;
      loop_ub = dIdx_data[b_i];
      if (b_d < loop_ub) {
        c_i = 0;
      } else {
        minSize = static_cast<int>(b_d) - loop_ub;
        c_i = minSize + 1;
        for (int i{0}; i <= minSize; i++) {
          y_data[i] = static_cast<double>(loop_ub) + static_cast<double>(i);
        }
      }
      minSize = idx[static_cast<int>(y_data[0]) - 1];
      for (int i{2}; i <= c_i; i++) {
        loop_ub = idx[static_cast<int>(y_data[i - 1]) - 1];
        if (minSize > loop_ub) {
          minSize = loop_ub;
        }
      }
      b_ib_data[b_i] = minSize;
    } else {
      loop_ub = dIdx_data[b_i];
      c_i = 7 - loop_ub;
      loop_ub = 6 - loop_ub;
      for (int i{0}; i <= loop_ub; i++) {
        ib_data[i] =
            static_cast<double>(dIdx_data[b_i]) + static_cast<double>(i);
      }
      minSize = idx[static_cast<int>(ib_data[0]) - 1];
      for (int i{2}; i <= c_i; i++) {
        loop_ub = idx[static_cast<int>(ib_data[i - 1]) - 1];
        if (minSize > loop_ub) {
          minSize = loop_ub;
        }
      }
      b_ib_data[b_i] = minSize;
    }
    loop_ub = anames[minSize - 1].f1.size[1];
    uB_data[b_i].f1.size[0] = 1;
    uB_data[b_i].f1.size[1] = loop_ub;
    for (int i{0}; i < loop_ub; i++) {
      uB_data[b_i].f1.data[i] = anames[b_ib_data[b_i] - 1].f1.data[i];
    }
  }
  std::memset(&uAInUB_data[0], 0, static_cast<unsigned int>(nz) * sizeof(bool));
  std::memset(&dIdx_data[0], 0, static_cast<unsigned int>(nz) * sizeof(int));
  j = 0;
  d_i = 0;
  exitg1 = false;
  while ((!exitg1) && (d_i <= nz - 1)) {
    int minSize_tmp;
    minSize_tmp = uA_data[d_i].f1.size[1];
    loop_ub = uB_data[j].f1.size[1];
    if (minSize_tmp <= loop_ub) {
      minSize = minSize_tmp;
    } else {
      minSize = loop_ub;
    }
    c_i = 0;
    while ((c_i + 1 <= minSize) &&
           (uA_data[d_i].f1.data[c_i] == uB_data[j].f1.data[c_i])) {
      c_i++;
    }
    if (c_i + 1 == minSize + 1) {
      is_less_than = (minSize_tmp < loop_ub);
    } else {
      is_less_than = (uA_data[d_i].f1.data[c_i] < uB_data[j].f1.data[c_i]);
    }
    if (!is_less_than) {
      if (internal::b_strcmp(uA_data[d_i].f1.data, uA_data[d_i].f1.size,
                             uB_data[j].f1.data, uB_data[j].f1.size)) {
        uAInUB_data[d_i] = true;
        dIdx_data[d_i] = b_ib_data[j];
        j++;
      } else {
        bool exitg2;
        exitg2 = false;
        while ((!exitg2) && (j + 1 <= b_nz)) {
          loop_ub = uB_data[j].f1.size[1];
          if (minSize_tmp <= loop_ub) {
            minSize = minSize_tmp;
          } else {
            minSize = loop_ub;
          }
          c_i = 0;
          while ((c_i + 1 <= minSize) &&
                 (uA_data[d_i].f1.data[c_i] == uB_data[j].f1.data[c_i])) {
            c_i++;
          }
          if (c_i + 1 == minSize + 1) {
            is_less_than = (minSize_tmp > loop_ub);
          } else {
            is_less_than =
                (uA_data[d_i].f1.data[c_i] > uB_data[j].f1.data[c_i]);
          }
          if (is_less_than) {
            j++;
          } else {
            exitg2 = true;
          }
        }
        if ((j + 1 <= b_nz) &&
            internal::b_strcmp(uA_data[d_i].f1.data, uA_data[d_i].f1.size,
                               uB_data[j].f1.data, uB_data[j].f1.size)) {
          uAInUB_data[d_i] = true;
          dIdx_data[d_i] = b_ib_data[j];
          j++;
        }
      }
    }
    if (j + 1 > b_nz) {
      exitg1 = true;
    } else {
      d_i++;
    }
  }
  for (int i{0}; i < 7; i++) {
    b2a[i] = 0U;
  }
  for (int i{0}; i < 6; i++) {
    loop_ub = static_cast<int>(icA[i]) - 1;
    minSize = dIdx_data[loop_ub];
    if (minSize < 0) {
      minSize = 0;
    } else if (minSize > 255) {
      minSize = 255;
    }
    b2a[i + 1] = static_cast<unsigned char>(minSize);
    d[i] = !uAInUB_data[loop_ub];
  }
  minSize = 0;
  loop_ub = 0;
  exitg1 = false;
  while ((!exitg1) && (loop_ub < 6)) {
    if (d[loop_ub]) {
      minSize++;
      ii_data[minSize - 1] = static_cast<signed char>(loop_ub + 1);
      if (minSize >= 6) {
        exitg1 = true;
      } else {
        loop_ub++;
      }
    } else {
      loop_ub++;
    }
  }
  if (minSize < 1) {
    minSize = 0;
  }
  for (int i{0}; i < minSize; i++) {
    ib_data[i] = ii_data[i];
  }
  is_less_than = true;
  loop_ub = 0;
  exitg1 = false;
  while ((!exitg1) && (loop_ub < 6)) {
    if (!uAInUB_data[static_cast<int>(icA[loop_ub]) - 1]) {
      is_less_than = false;
      exitg1 = true;
    } else {
      loop_ub++;
    }
  }
  if (!is_less_than) {
    if (minSize < 1) {
      y.set_size(1, 0);
    } else {
      y.set_size(1, minSize);
      for (int i{0}; i < minSize; i++) {
        y[i] = static_cast<signed char>(i + 1);
      }
    }
    loop_ub = y.size(1);
    for (int i{0}; i < loop_ub; i++) {
      b2a[static_cast<int>(ib_data[i])] =
          static_cast<unsigned char>(static_cast<unsigned int>(y[i]) + 6U);
    }
  }
  for (int i{0}; i < 6; i++) {
    bcodesout[i] = b2a[i + 1];
  }
}

} // namespace coder

// End of code generation (categorical.cpp)
