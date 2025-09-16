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
#include "rt_nonfinite.h"
#include "sortIdx.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void sort(double x_data[], const int *x_size)
{
  double vwork_data[150];
  int dim;
  int j;
  int k;
  int vstride;
  int vwork_size;
  dim = 2;
  if (*x_size != 1) {
    dim = 1;
  }
  if (dim <= 1) {
    vwork_size = *x_size;
  } else {
    vwork_size = 1;
  }
  vstride = 1;
  dim -= 2;
  for (k = 0; k <= dim; k++) {
    vstride *= *x_size;
  }
  for (j = 0; j < vstride; j++) {
    double xwork_data[150];
    double x4[4];
    int iidx_data[150];
    int iwork_data[150];
    int i1;
    int i2;
    int i4;
    int ib;
    int nNaNs;
    unsigned char idx4[4];
    for (k = 0; k < vwork_size; k++) {
      vwork_data[k] = x_data[j + k * vstride];
      iidx_data[k] = 0;
      iwork_data[k] = 0;
    }
    x4[0] = 0.0;
    idx4[0] = 0U;
    x4[1] = 0.0;
    idx4[1] = 0U;
    x4[2] = 0.0;
    idx4[2] = 0U;
    x4[3] = 0.0;
    idx4[3] = 0U;
    nNaNs = 0;
    ib = 0;
    for (k = 0; k < vwork_size; k++) {
      if (rtIsNaN(vwork_data[k])) {
        dim = (vwork_size - nNaNs) - 1;
        iidx_data[dim] = k + 1;
        xwork_data[dim] = vwork_data[k];
        nNaNs++;
      } else {
        ib++;
        idx4[ib - 1] = (unsigned char)(k + 1);
        x4[ib - 1] = vwork_data[k];
        if (ib == 4) {
          double d;
          double d1;
          int b_i1;
          int i;
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
              i = i1;
              b_i1 = i2;
              i1 = ib;
              i2 = i4;
            } else if (x4[i2 - 1] <= x4[i4 - 1]) {
              i = i1;
              b_i1 = ib;
              i1 = i2;
              i2 = i4;
            } else {
              i = i1;
              b_i1 = ib;
              i1 = i4;
            }
          } else if (d <= x4[i4 - 1]) {
            if (x4[i2 - 1] <= x4[i4 - 1]) {
              i = ib;
              b_i1 = i1;
              i1 = i2;
              i2 = i4;
            } else {
              i = ib;
              b_i1 = i1;
              i1 = i4;
            }
          } else {
            i = ib;
            b_i1 = i4;
          }
          iidx_data[dim - 3] = idx4[i - 1];
          iidx_data[dim - 2] = idx4[b_i1 - 1];
          iidx_data[dim - 1] = idx4[i1 - 1];
          iidx_data[dim] = idx4[i2 - 1];
          vwork_data[dim - 3] = x4[i - 1];
          vwork_data[dim - 2] = x4[b_i1 - 1];
          vwork_data[dim - 1] = x4[i1 - 1];
          vwork_data[dim] = x4[i2 - 1];
          ib = 0;
        }
      }
    }
    i4 = vwork_size - nNaNs;
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
        i1 = (i4 - ib) + k;
        i2 = perm[k];
        iidx_data[i1] = idx4[i2 - 1];
        vwork_data[i1] = x4[i2 - 1];
      }
    }
    dim = nNaNs >> 1;
    for (k = 0; k < dim; k++) {
      i1 = i4 + k;
      i2 = iidx_data[i1];
      ib = (vwork_size - k) - 1;
      iidx_data[i1] = iidx_data[ib];
      iidx_data[ib] = i2;
      vwork_data[i1] = xwork_data[ib];
      vwork_data[ib] = xwork_data[i1];
    }
    if (((unsigned int)nNaNs & 1U) != 0U) {
      dim += i4;
      vwork_data[dim] = xwork_data[dim];
    }
    if (i4 > 1) {
      i2 = i4 >> 2;
      ib = 4;
      while (i2 > 1) {
        if (((unsigned int)i2 & 1U) != 0U) {
          i2--;
          dim = ib * i2;
          i1 = i4 - dim;
          if (i1 > ib) {
            merge(iidx_data, vwork_data, dim, ib, i1 - ib, iwork_data,
                  xwork_data);
          }
        }
        dim = ib << 1;
        i2 >>= 1;
        for (k = 0; k < i2; k++) {
          merge(iidx_data, vwork_data, k * dim, ib, ib, iwork_data, xwork_data);
        }
        ib = dim;
      }
      if (i4 > ib) {
        merge(iidx_data, vwork_data, 0, ib, i4 - ib, iwork_data, xwork_data);
      }
    }
    for (k = 0; k < vwork_size; k++) {
      x_data[j + k * vstride] = vwork_data[k];
    }
  }
}

/* End of code generation (sort.c) */
