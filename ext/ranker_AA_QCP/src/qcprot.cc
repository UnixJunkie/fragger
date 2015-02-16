/******************************************************************************
 *
 * Function: Rapid calculation of the least-squares rotation using a
 * quaternion-based characteristic polynomial and a cofactor matrix
 *
 * Author(s): Douglas L. Theobald
 * Department of Biochemistry
 * MS 009
 * Brandeis University
 * 415 South St
 * Waltham, MA 02453
 * USA
 *
 * dtheobald@brandeis.edu
 *
 * Pu Liu
 * Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
 * 665 Stockton Drive
 * Exton, PA 19341
 * USA
 *
 * pliu24@its.jnj.com
 *
 * Tweaks for use in Durandal by Francois Berenger, RIKEN, Japan.
 *
 * If you use this QCP rotation calculation method in a publication, please
 * reference:
 *
 * Douglas L. Theobald (2005)
 * "Rapid calculation of RMSD using a quaternion-based characteristic
 * polynomial." Acta Crystallographica A 61(4):478-480.
 *
 * Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
 * "Fast determination of the optimal rotational matrix for macromolecular
 * superpositions." in press, Journal of Computational Chemistry
 *
 * Copyright (c) 2009-2010, Pu Liu and Douglas L. Theobald All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of the <ORGANIZATION> nor the names of its contributors
 *   may be used to endorse or promote products derived from this software
 *   without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *****************************************************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <limits>

#include "qcprot.h"

double **MatInit(const int rows, const int cols) {

  int             i;
  double        **matrix = NULL;
  double         *matspace = NULL;

  matspace = (double *) calloc((rows * cols), sizeof(double));
  if (matspace == NULL) {
    perror("\n ERROR");
    printf("\n ERROR: Failure to allocate matrix space in "
           "MatInit(): (%d x %d)\n", rows, cols);
    exit(EXIT_FAILURE);
  }

  /* allocate room for the pointers to the rows */
  matrix = (double **) malloc(rows * sizeof(double *));
  if (matrix == NULL) {
    perror("\n ERROR");
    printf("\n ERROR: Failure to allocate room for row pointers in "
           "MatInit(): (%d)\n", rows);
    exit(EXIT_FAILURE);
  }

  /*  now 'point' the pointers */
  for (i = 0; i < rows; i++)
    matrix[i] = matspace + (i * cols);

  return(matrix);
}

void MatDestroy(double ***matrix_ptr) {

  double **matrix = *matrix_ptr;

  if (matrix != NULL) {
    if (matrix[0] != NULL) {
      free(matrix[0]);
      matrix[0] = NULL;
    }
    free(matrix);
    *matrix_ptr = NULL;
  }
}

double G1 = 0.0;

double InnerProduct(double *A,
                    double **coords1,
                    double **coords2,
                    const int len,
                    bool update_G1) {

  double          x1, x2, y1, y2, z1, z2;
  int             i;
  const double   *fx1 = coords1[0], *fy1 = coords1[1], *fz1 = coords1[2];
  const double   *fx2 = coords2[0], *fy2 = coords2[1], *fz2 = coords2[2];
  double         G2;

  G2 = A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

  if (update_G1) {
    G1 = 0.0;
    for (i = 0; i < len; ++i) {

      x1 = fx1[i];
      y1 = fy1[i];
      z1 = fz1[i];
      x2 = fx2[i];
      y2 = fy2[i];
      z2 = fz2[i];

      G1 += x1 * x1 + y1 * y1 + z1 * z1;
      G2 += x2 * x2 + y2 * y2 + z2 * z2;

      A[0] +=  (x1 * x2);
      A[1] +=  (x1 * y2);
      A[2] +=  (x1 * z2);
      A[3] +=  (y1 * x2);
      A[4] +=  (y1 * y2);
      A[5] +=  (y1 * z2);
      A[6] +=  (z1 * x2);
      A[7] +=  (z1 * y2);
      A[8] +=  (z1 * z2);
    }
  } else {
    for (i = 0; i < len; ++i) {

      x1 = fx1[i];
      y1 = fy1[i];
      z1 = fz1[i];
      x2 = fx2[i];
      y2 = fy2[i];
      z2 = fz2[i];

      G2 += x2 * x2 + y2 * y2 + z2 * z2;

      A[0] +=  (x1 * x2);
      A[1] +=  (x1 * y2);
      A[2] +=  (x1 * z2);
      A[3] +=  (y1 * x2);
      A[4] +=  (y1 * y2);
      A[5] +=  (y1 * z2);
      A[6] +=  (z1 * x2);
      A[7] +=  (z1 * y2);
      A[8] +=  (z1 * z2);
    }
  }

  return (G1 + G2) * 0.5;
}

double FastCalcRMSD(double *A,
                    double E0,
                    int len,
                    size_t idx1,
                    size_t idx2) {

  double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
  double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
    SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
    SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
    SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
  double C[4];
  int i;
  double mxEigenV;
  double oldg = 0.0;
  double b, a, delta, x2;
  double evalprec = 1e-9;

  Sxx = A[0]; Sxy = A[1]; Sxz = A[2];
  Syx = A[3]; Syy = A[4]; Syz = A[5];
  Szx = A[6]; Szy = A[7]; Szz = A[8];

  Sxx2 = Sxx * Sxx;
  Syy2 = Syy * Syy;
  Szz2 = Szz * Szz;

  Sxy2 = Sxy * Sxy;
  Syz2 = Syz * Syz;
  Sxz2 = Sxz * Sxz;

  Syx2 = Syx * Syx;
  Szy2 = Szy * Szy;
  Szx2 = Szx * Szx;

  SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
  Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

  C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
  C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx -
                Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

  SxzpSzx = Sxz + Szx;
  SyzpSzy = Syz + Szy;
  SxypSyx = Sxy + Syx;
  SyzmSzy = Syz - Szy;
  SxzmSzx = Sxz - Szx;
  SxymSyx = Sxy - Syx;
  SxxpSyy = Sxx + Syy;
  SxxmSyy = Sxx - Syy;
  Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

  C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
    + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) *
    (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
    + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) *
    (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
    + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) *
    (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
    + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) *
    (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
    + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) *
    (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

  mxEigenV = E0;
  int max_trials = 50;
  for (i = 0; i < max_trials; ++i) {
    oldg = mxEigenV;
    x2 = mxEigenV*mxEigenV;
    b = (x2 + C[2])*mxEigenV;
    a = b + C[1];
    delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
    mxEigenV -= delta;
    if (fabs(mxEigenV - oldg) < fabs((evalprec)*mxEigenV))
      break;
  }

  if (i == max_trials) {
    fprintf(stderr, "%s:%d:%s: did not converge for %ld %ld\n",
            __FILE__, __LINE__, __FUNCTION__, idx1, idx2);
  }

  double E0_minus_mxEigenV = E0 - mxEigenV;
  if (E0_minus_mxEigenV < 0.0) {
    // avoid NaN case when RMSD of 2 identical structures is being computed
    // and mxEigenV is slightly larger than E0
    return 0.0;
  } else {
    return sqrt(2.0 * E0_minus_mxEigenV / len);
  }
}

double rmsd_without_rotation_matrix(double **coords1,
                                    double **coords2,
                                    const int len,
                                    bool different_ref,
                                    size_t i, size_t j) {

  double A[9];
  /* center the structures */
  //CenterCoords(coords1, len); // in Durandal PDBs are centered
  //CenterCoords(coords2, len); // once and for all at read time

  /* calculate the inner product of two structures */
  double E0 = InnerProduct(A, coords1, coords2, len, different_ref);

  return FastCalcRMSD(A, E0, len, i, j);
}
