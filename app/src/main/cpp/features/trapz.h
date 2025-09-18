/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * trapz.h
 *
 * Code generation for function 'trapz'
 *
 */

#ifndef TRAPZ_H
#define TRAPZ_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
double trapz(const double x_data[], int x_size, const double y_data[],
             int y_size);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (trapz.h) */
