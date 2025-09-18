/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CompactClassificationModel.h
 *
 * Code generation for function 'CompactClassificationModel'
 *
 */

#ifndef COMPACTCLASSIFICATIONMODEL_H
#define COMPACTCLASSIFICATIONMODEL_H

/* Include files */
#include "predict_exercise_internal_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void c_CompactClassificationModel_Co(c_classreg_learning_classif_Com *obj);

int c_CompactClassificationModel_mi(const double obj_Cost[36],
                                    const double scoresIn_data[],
                                    const int scoresIn_size[2],
                                    double classnum_data[], double cost_data[],
                                    int cost_size[2]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (CompactClassificationModel.h) */
