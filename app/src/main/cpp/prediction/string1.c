/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * string1.c
 *
 * Code generation for function 'string1'
 *
 */

/* Include files */
#include "string1.h"
#include "rt_nonfinite.h"

/* Function Definitions */
boolean_T string_gt(const char obj1_Value_data[], const int obj1_Value_size[2],
                    const char obj2_Value_data[], const int obj2_Value_size[2])
{
  int i;
  int minSize;
  boolean_T is_greater_than;
  i = obj1_Value_size[1];
  minSize = obj2_Value_size[1];
  if (i <= minSize) {
    minSize = i;
  }
  if (minSize == 0) {
    is_greater_than = (obj1_Value_size[1] > obj2_Value_size[1]);
  } else {
    i = 0;
    while ((i + 1 <= minSize) && (obj1_Value_data[i] == obj2_Value_data[i])) {
      i++;
    }
    if (i + 1 == minSize + 1) {
      is_greater_than = (obj1_Value_size[1] > obj2_Value_size[1]);
    } else {
      is_greater_than = (obj1_Value_data[i] > obj2_Value_data[i]);
    }
  }
  return is_greater_than;
}

boolean_T string_lt(const char obj1_Value_data[], const int obj1_Value_size[2],
                    const char obj2_Value_data[], const int obj2_Value_size[2])
{
  int i;
  int minSize;
  boolean_T is_less_than;
  i = obj1_Value_size[1];
  minSize = obj2_Value_size[1];
  if (i <= minSize) {
    minSize = i;
  }
  if (minSize == 0) {
    is_less_than = (obj1_Value_size[1] < obj2_Value_size[1]);
  } else {
    i = 0;
    while ((i + 1 <= minSize) && (obj1_Value_data[i] == obj2_Value_data[i])) {
      i++;
    }
    if (i + 1 == minSize + 1) {
      is_less_than = (obj1_Value_size[1] < obj2_Value_size[1]);
    } else {
      is_less_than = (obj1_Value_data[i] < obj2_Value_data[i]);
    }
  }
  return is_less_than;
}

/* End of code generation (string1.c) */
