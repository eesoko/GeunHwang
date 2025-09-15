//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// periodogram.h
//
// Code generation for function 'periodogram'
//

#ifndef PERIODOGRAM_H
#define PERIODOGRAM_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void periodogram(const array<double, 1U> &x, array<double, 1U> &Px,
                 array<double, 1U> &w);

}

#endif
// End of code generation (periodogram.h)
