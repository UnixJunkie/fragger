/* **************************************************************************
 * Copyright 2009 Shuai Cheng Li and Yen Kaow Ng
 * **************************************************************************
 * In 2012 and 2013 Francois Berenger made many changes
 * (more robust, new constructors, binary mode, etc.).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * **************************************************************************/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;

#include "SimpPDB.h"

SimPDB::SimPDB(vector<string>& fragment_lines,
               vector<double>& coordinates,
               int expected_count,
               double** dst) {

  mNumResidue = 0; // we are in all atoms mode so this is an atom count in fact
  int count = 0;

  mProteinFileName = fragment_lines[0];
  if (mProteinFileName.substr(0, 3) != "TER") {
    // we expect the first line of a fragment to be TER followed by fragment id
    cerr << __FILE__ << ": " << __LINE__
         << ": ERROR: non TER line starting fragment in: "
         << mProteinFileName << endl;
    exit(1);
  }

  if (fragment_lines.size() > 1) { // text mode

    for (size_t i = 1; i < fragment_lines.size(); ++i) {

      if (fragment_lines[i].substr(0, 4) != "ATOM") {
        cerr << __FILE__ << ": " << __LINE__
             << ": ERROR: non ATOM line in fragment: "
             << fragment_lines[i] << endl;
        exit(1);
      } else {
        // next line not faster than 3 atof calls
        //sscanf(fragment_lines[i].substr(30, 24).c_str(),"%f%f%f",&x, &y, &z);
        dst[0][count] = atof(fragment_lines[i].substr(30, 8).c_str()); // x
        dst[1][count] = atof(fragment_lines[i].substr(38, 8).c_str()); // y
        dst[2][count] = atof(fragment_lines[i].substr(46, 8).c_str()); // z
        ++count;
      }
    }

  } else if (coordinates.size() > 0) { // binary mode

    if (coordinates.size() % 3 != 0) {
      cerr << __FILE__ << ": " << __LINE__
           << ": ERROR: size of coordinates: "
           << coordinates.size() << " in " << mProteinFileName
           << endl;
      exit(1);
    }

    for (count = 0; count < ((int)coordinates.size()) / 3; ++count) {
      int i = 3 * count;
      dst[0][count] = coordinates[i    ]; // x
      dst[1][count] = coordinates[i + 1]; // y
      dst[2][count] = coordinates[i + 2]; // z
    }

  } else {
    cerr << __FILE__ << ": " << __LINE__
         << ": ERROR: no fragment_lines and no coordinates for "
         << mProteinFileName << endl;
    exit(1);
  }

  mNumResidue = count;
  if (count != expected_count) {
    cerr << __FILE__ << ": " << __LINE__
         << ": FATAL: " << mProteinFileName << ": read " << count
         << " atoms while expecting " << expected_count << endl;
    exit(1);
  }

  // centering (center of mass will be at (0, 0, 0) after)
  float cx, cy, cz;
  cx = cy = cz = 0;

  for(int i = 0; i < mNumResidue; ++i) {
    cx += dst[0][i];
    cy += dst[1][i];
    cz += dst[2][i];
  }

  cx /= mNumResidue;
  cy /= mNumResidue;
  cz /= mNumResidue;

  for(int i = 0; i < mNumResidue; ++i) {
    dst[0][i] -= cx;
    dst[1][i] -= cy;
    dst[2][i] -= cz;
  }
}
