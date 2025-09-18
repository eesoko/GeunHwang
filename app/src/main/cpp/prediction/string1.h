/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * string1.h
 *
 * Code generation for function 'string1'
 *
 */

#ifndef STRING1_H
#define STRING1_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
boolean_T string_gt(const char obj1_Value_data[], const int obj1_Value_size[2],
                    const char obj2_Value_data[], const int obj2_Value_size[2]);

boolean_T string_lt(const char obj1_Value_data[], const int obj1_Value_size[2],
                    const char obj2_Value_data[], const int obj2_Value_size[2]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (string1.h) */
