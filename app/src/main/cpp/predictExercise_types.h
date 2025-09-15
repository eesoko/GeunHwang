//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// predictExercise_types.h
//
// Code generation for function 'predictExercise'
//

#ifndef PREDICTEXERCISE_TYPES_H
#define PREDICTEXERCISE_TYPES_H

// Include files
#include "rtwtypes.h"
#include "coder_bounded_array.h"
#define MAX_THREADS omp_get_max_threads()

// Type Definitions
struct cell_wrap_3 {
  coder::bounded_array<char, 18U, 2U> f1;
};

#endif
// End of code generation (predictExercise_types.h)
