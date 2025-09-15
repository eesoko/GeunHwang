//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// predictExercise_terminate.cpp
//
// Code generation for function 'predictExercise_terminate'
//

// Include files
#include "predictExercise_terminate.h"
#include "predictExercise_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Definitions
void predictExercise_terminate()
{
  omp_destroy_nest_lock(&predictExercise_nestLockGlobal);
  isInitialized_predictExercise = false;
}

// End of code generation (predictExercise_terminate.cpp)
