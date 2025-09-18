/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CholeskyDiscriminant.h
 *
 * Code generation for function 'CholeskyDiscriminant'
 *
 */

#ifndef CHOLESKYDISCRIMINANT_H
#define CHOLESKYDISCRIMINANT_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void CholeskyDiscriminant(double mah_data[], const int mah_size[2],
                          const boolean_T classLogicalIndices[6],
                          const double X_data[], const int X_size[2],
                          const double Mu[96], const double invD[16],
                          const double invR[256]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (CholeskyDiscriminant.h) */
