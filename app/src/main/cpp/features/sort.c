/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sort.c
 *
 * Code generation for function 'sort'
 *
 */

/* Include files */
#include "sort.h"
#include "feature_extractor_codegen_emxutil.h"
#include "feature_extractor_codegen_types.h"
#include "rt_nonfinite.h"
#include "sortIdx.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void sort(emxArray_real_T *x)
{
  emxArray_int32_T *iidx;
  emxArray_int32_T *iwork;
  emxArray_real_T *vwork;
  double *vwork_data;
  double *x_data;
  double *xwork_data;
  int b;
  int b_b;
  int b_k;
  int dim;
  int i;
  int i1;
  int j;
  int k;
  int vstride;
  int *iidx_data;
  int *iwork_data;
  x_data = x->data;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }
  if (dim <= 1) {
    i = x->size[0];
  } else {
    i = 1;
  }
  emxInit_real_T(&vwork, 1);
  i1 = vwork->size[0];
  vwork->size[0] = i;
  emxEnsureCapacity_real_T(vwork, i1);
  vwork_data = vwork->data;
  vstride = 1;
  i1 = dim - 2;
  for (k = 0; k <= i1; k++) {
    vstride *= x->size[0];
  }
  emxInit_int32_T(&iidx, 1);
  emxInit_int32_T(&iwork, 1);
  emxInit_real_T(&x, 1);
  for (j = 0; j < vstride; j++) {
    double x4[4];
    int idx4[4];
    int bLen;
    int bLen2;
    int i2;
    int i4;
    int ib;
    int nBlocks;
    int nNaNs;
    int wOffset;
    for (k = 0; k < i; k++) {
      vwork_data[k] = x_data[j + k * vstride];
    }
    nBlocks = vwork->size[0];
    i1 = iidx->size[0];
    iidx->size[0] = vwork->size[0];
    emxEnsureCapacity_int32_T(iidx, i1);
    iidx_data = iidx->data;
    for (k = 0; k < nBlocks; k++) {
      iidx_data[k] = 0;
    }
    i1 = iwork->size[0];
    iwork->size[0] = vwork->size[0];
    emxEnsureCapacity_int32_T(iwork, i1);
    iwork_data = iwork->data;
    for (k = 0; k < nBlocks; k++) {
      iwork_data[k] = 0;
    }
    i1 = x->size[0];
    x->size[0] = vwork->size[0];
    emxEnsureCapacity_real_T(x, i1);
    xwork_data = x->data;
    x4[0] = 0.0;
    idx4[0] = 0;
    x4[1] = 0.0;
    idx4[1] = 0;
    x4[2] = 0.0;
    idx4[2] = 0;
    x4[3] = 0.0;
    idx4[3] = 0;
    nNaNs = 0;
    ib = 0;
    for (k = 0; k < nBlocks; k++) {
      if (rtIsNaN(vwork_data[k])) {
        dim = (nBlocks - nNaNs) - 1;
        iidx_data[dim] = k + 1;
        xwork_data[dim] = vwork_data[k];
        nNaNs++;
      } else {
        ib++;
        idx4[ib - 1] = k + 1;
        x4[ib - 1] = vwork_data[k];
        if (ib == 4) {
          double d;
          double d1;
          dim = k - nNaNs;
          if (x4[0] <= x4[1]) {
            i1 = 1;
            i2 = 2;
          } else {
            i1 = 2;
            i2 = 1;
          }
          if (x4[2] <= x4[3]) {
            ib = 3;
            i4 = 4;
          } else {
            ib = 4;
            i4 = 3;
          }
          d = x4[i1 - 1];
          d1 = x4[ib - 1];
          if (d <= d1) {
            if (x4[i2 - 1] <= d1) {
              bLen = i1;
              bLen2 = i2;
              i1 = ib;
              i2 = i4;
            } else if (x4[i2 - 1] <= x4[i4 - 1]) {
              bLen = i1;
              bLen2 = ib;
              i1 = i2;
              i2 = i4;
            } else {
              bLen = i1;
              bLen2 = ib;
              i1 = i4;
            }
          } else if (d <= x4[i4 - 1]) {
            if (x4[i2 - 1] <= x4[i4 - 1]) {
              bLen = ib;
              bLen2 = i1;
              i1 = i2;
              i2 = i4;
            } else {
              bLen = ib;
              bLen2 = i1;
              i1 = i4;
            }
          } else {
            bLen = ib;
            bLen2 = i4;
          }
          iidx_data[dim - 3] = idx4[bLen - 1];
          iidx_data[dim - 2] = idx4[bLen2 - 1];
          iidx_data[dim - 1] = idx4[i1 - 1];
          iidx_data[dim] = idx4[i2 - 1];
          vwork_data[dim - 3] = x4[bLen - 1];
          vwork_data[dim - 2] = x4[bLen2 - 1];
          vwork_data[dim - 1] = x4[i1 - 1];
          vwork_data[dim] = x4[i2 - 1];
          ib = 0;
        }
      }
    }
    wOffset = vwork->size[0] - nNaNs;
    if (ib > 0) {
      signed char perm[4];
      perm[1] = 0;
      perm[2] = 0;
      perm[3] = 0;
      if (ib == 1) {
        perm[0] = 1;
      } else if (ib == 2) {
        if (x4[0] <= x4[1]) {
          perm[0] = 1;
          perm[1] = 2;
        } else {
          perm[0] = 2;
          perm[1] = 1;
        }
      } else if (x4[0] <= x4[1]) {
        if (x4[1] <= x4[2]) {
          perm[0] = 1;
          perm[1] = 2;
          perm[2] = 3;
        } else if (x4[0] <= x4[2]) {
          perm[0] = 1;
          perm[1] = 3;
          perm[2] = 2;
        } else {
          perm[0] = 3;
          perm[1] = 1;
          perm[2] = 2;
        }
      } else if (x4[0] <= x4[2]) {
        perm[0] = 2;
        perm[1] = 1;
        perm[2] = 3;
      } else if (x4[1] <= x4[2]) {
        perm[0] = 2;
        perm[1] = 3;
        perm[2] = 1;
      } else {
        perm[0] = 3;
        perm[1] = 2;
        perm[2] = 1;
      }
      dim = (unsigned char)ib;
      for (k = 0; k < dim; k++) {
        i1 = (wOffset - ib) + k;
        i2 = perm[k];
        iidx_data[i1] = idx4[i2 - 1];
        vwork_data[i1] = x4[i2 - 1];
      }
    }
    dim = nNaNs >> 1;
    for (k = 0; k < dim; k++) {
      i1 = wOffset + k;
      i2 = iidx_data[i1];
      ib = (nBlocks - k) - 1;
      iidx_data[i1] = iidx_data[ib];
      iidx_data[ib] = i2;
      vwork_data[i1] = xwork_data[ib];
      vwork_data[ib] = xwork_data[i1];
    }
    if (((unsigned int)nNaNs & 1U) != 0U) {
      i1 = wOffset + dim;
      vwork_data[i1] = xwork_data[i1];
    }
    dim = 2;
    if (wOffset > 1) {
      if (vwork->size[0] >= 256) {
        nBlocks = wOffset >> 8;
        if (nBlocks > 0) {
          for (b = 0; b < nBlocks; b++) {
            double xwork[256];
            int b_iwork[256];
            i4 = (b << 8) - 1;
            for (b_b = 0; b_b < 6; b_b++) {
              bLen = 1 << (b_b + 2);
              bLen2 = bLen << 1;
              nNaNs = 256 >> (b_b + 3);
              for (b_k = 0; b_k < nNaNs; b_k++) {
                i1 = (i4 + b_k * bLen2) + 1;
                for (k = 0; k < bLen2; k++) {
                  dim = i1 + k;
                  b_iwork[k] = iidx_data[dim];
                  xwork[k] = vwork_data[dim];
                }
                i2 = 0;
                ib = bLen;
                dim = i1 - 1;
                int exitg1;
                do {
                  exitg1 = 0;
                  dim++;
                  if (xwork[i2] <= xwork[ib]) {
                    iidx_data[dim] = b_iwork[i2];
                    vwork_data[dim] = xwork[i2];
                    if (i2 + 1 < bLen) {
                      i2++;
                    } else {
                      exitg1 = 1;
                    }
                  } else {
                    iidx_data[dim] = b_iwork[ib];
                    vwork_data[dim] = xwork[ib];
                    if (ib + 1 < bLen2) {
                      ib++;
                    } else {
                      dim -= i2;
                      for (k = i2 + 1; k <= bLen; k++) {
                        ib = dim + k;
                        iidx_data[ib] = b_iwork[k - 1];
                        vwork_data[ib] = xwork[k - 1];
                      }
                      exitg1 = 1;
                    }
                  }
                } while (exitg1 == 0);
              }
            }
          }
          dim = nBlocks << 8;
          i1 = wOffset - dim;
          if (i1 > 0) {
            merge_block(iidx, vwork, dim, i1, 2, iwork, x);
          }
          dim = 8;
        }
      }
      merge_block(iidx, vwork, 0, wOffset, dim, iwork, x);
      vwork_data = vwork->data;
    }
    for (k = 0; k < i; k++) {
      x_data[j + k * vstride] = vwork_data[k];
    }
  }
  emxFree_real_T(&x);
  emxFree_int32_T(&iwork);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

/* End of code generation (sort.c) */
