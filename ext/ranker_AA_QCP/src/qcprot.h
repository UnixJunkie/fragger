/**
 *  Copyright (C) 2013, Zhang Initiative Research Unit,
 *  Advance Science Institute, Riken
 *  2-1 Hirosawa, Wako, Saitama 351-0198, Japan
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************************/

#ifndef QCPROT_H
#define QCPROT_H

double **MatInit(const int rows, const int cols);

void MatDestroy(double ***matrix_ptr);

double rmsd_without_rotation_matrix(double **coords1,
                                    double **coords2,
                                    const int len,
                                    bool different_ref,
                                    size_t i,
                                    size_t j);

#endif /* QCPROT_H */
