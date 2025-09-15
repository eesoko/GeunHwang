//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// categorical.h
//
// Code generation for function 'categorical'
//

#ifndef CATEGORICAL_H
#define CATEGORICAL_H

// Include files
#include "predictExercise_types.h"
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Definitions
namespace coder {
class categorical {
public:
  static void convertCodes(const cell_wrap_3 bnames[6],
                           const cell_wrap_3 anames[6],
                           unsigned char bcodesout[6]);
  categorical();
  ~categorical();
  unsigned char codes;
  cell_wrap_3 categoryNames[6];
};

} // namespace coder

#endif
// End of code generation (categorical.h)
