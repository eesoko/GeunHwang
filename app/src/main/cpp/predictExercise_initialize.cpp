//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// predictExercise_initialize.cpp
//
// Code generation for function 'predictExercise_initialize'
//

// Include files
#include "predictExercise_initialize.h"
#include "predictExercise.h"
#include "predictExercise_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Definitions
void predictExercise_initialize()
{
  omp_init_nest_lock(&predictExercise_nestLockGlobal);
  predictExercise_init();
  isInitialized_predictExercise = true;
}

// End of code generation (predictExercise_initialize.cpp)
