/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * cellstr_unique.h
 *
 * Code generation for function 'cellstr_unique'
 *
 */

#ifndef CELLSTR_UNIQUE_H
#define CELLSTR_UNIQUE_H

/* Include files */
#include "predict_exercise_internal_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
int b_cellstr_unique(const cell_wrap_3 a[6], cell_wrap_3 u_data[]);

int cellstr_unique(const cell_wrap_3 a[6], cell_wrap_3 u_data[],
                   double ia_data[], int *ia_size);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (cellstr_unique.h) */
