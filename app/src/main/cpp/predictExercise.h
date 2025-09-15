//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// predictExercise.h
//
// Code generation for function 'predictExercise'
//

#ifndef PREDICTEXERCISE_H
#define PREDICTEXERCISE_H

// Include files
#include "categorical.h"
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void predictExercise(const coder::array<double, 2U> &sensorData,
                            coder::categorical *predictedExercise);

void predictExercise_init();

#endif
// End of code generation (predictExercise.h)
