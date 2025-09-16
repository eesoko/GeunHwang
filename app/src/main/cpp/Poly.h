/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Poly.h
 *
 * Code generation for function 'Poly'
 *
 */

#ifndef POLY_H
#define POLY_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void Poly(const double svT_data[], const int svT_size[2], const double x[32],
          double kernelProduct_data[], int kernelProduct_size[2]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (Poly.h) */
