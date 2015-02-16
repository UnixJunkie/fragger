/**
 * Copyright (C) 2013, Zhang Initiative Research Unit,
 * Advance Science Institute, Riken
 * 2-1 Hirosawa, Wako, Saitama 351-0198, Japan
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
 *****************************************************************************/

#include <arpa/inet.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "SimpPDB.h"
#include "qcprot.h"

int _mLen        = 0;
size_t _previous = (size_t) -1;
double** _mPOS1  = NULL;
double** _mPOS2  = NULL;

// WARNING: this function is not thread safe, not suitable for OpenMP
float rmsd(size_t idx1, size_t idx2) {

  float rmsd = rmsd_without_rotation_matrix(_mPOS1, _mPOS2, _mLen,
                                            idx1 != _previous,
                                            idx1, idx2);
  _previous = idx1;

  return rmsd;
}

float lower_bound_rmsd() {

  float x1, y1, z1, x2, y2, z2, nui, nvi, acc = 0.0;

  for (int i = 0 ; i < _mLen ; ++i) {
    x1 = fabs(_mPOS1[0][i]);
    y1 = fabs(_mPOS1[1][i]);
    z1 = fabs(_mPOS1[2][i]);
    x1 *= x1;
    y1 *= y1;
    z1 *= z1;
    nui = sqrt(x1 + y1 + z1);

    x2 = fabs(_mPOS2[0][i]);
    y2 = fabs(_mPOS2[1][i]);
    z2 = fabs(_mPOS2[2][i]);
    x2 *= x2;
    y2 *= y2;
    z2 *= z2;
    nvi = sqrt(x2 + y2 + z2);

    nui -= nvi;
    nui *= nui;
    acc += nui;
  }

  return sqrt(acc / _mLen);
}

// transform "1,2,3" into <1,2,3>
void extract_res_nums(string& s, vector<int>& out) {

  string delimiters = ",";
  size_t current;
  size_t next = -1;

  out.clear();
  do {
    current = next + 1;
    next = s.find_first_of( delimiters, current );
    out.push_back( atoi(s.substr( current, next - current ).c_str()) );
  } while (next != string::npos);
}

// keep first line followed by lines of the selected residues only
void reduce(vector<string>& fragment_lines_in,
            vector<string>& fragment_lines_out,
            vector<int>& selected_residues) {

  fragment_lines_out.clear();
  fragment_lines_out.push_back(fragment_lines_in[0]); // copy TER line

  for (size_t i = 0; i < selected_residues.size(); ++i) {

    int res_num = selected_residues[i];
    int j = 1 + ((res_num - 1) * 4);
    // copy the 4 bb atoms
    fragment_lines_out.push_back(fragment_lines_in[j + 0]);
    fragment_lines_out.push_back(fragment_lines_in[j + 1]);
    fragment_lines_out.push_back(fragment_lines_in[j + 2]);
    fragment_lines_out.push_back(fragment_lines_in[j + 3]);
  }
}

// read (no decoding) the fragment header + its binary coordinates
// return true if more fragment remain to be read after
bool raw_read_one_binary_fragment(int read_size,
                                  ifstream& input,
                                  string& header,
                                  vector<unsigned char>& payload) {
  header.clear();
  payload.clear();

  header.reserve(80);
  payload.reserve(read_size);
  // read header
  getline(input, header);
  // read payload
  unsigned char c;
  for (int i = 0 ; i < read_size ; ++i) {
    c = input.get();
    payload.push_back(c);
  }

  return (input.peek() != EOF);
}

int main (int argc, char** argv) {

  vector<string> ref_frag;
  vector<double> empty_coords;
  vector<int> selected_residues;

  if (argc < 5 || argc > 6) {
    cout << "Output on stdout the All Atom RMSD of each PDB fragment\n"
         << "(in the fragments file) to the reference fragment\n"
         << "if -bin is given at the end, fragments are read\n"
         << "in binary format\n"
         << "rmsd_limit (float): if not 0: don't output fragments above "
         << "this RMSD\n"
         << "---" << endl;
    //              0        1                2              3
    cout << "usage: ./ranker {nb_res|1,2,7,8} fragments_file ref_frag "
      //     4          5
         << "rmsd_limit [-bin]"
         << endl;
    return 1;
  }

  // if there is at least one quote in argv[1]: use extract_res_nums
  int n = 0;
  int m = -1; // -1 because we don't want to count the TER line
  string argv_1_str = argv[1];
  string delimiters = ",";
  size_t current = 0;
  if (argv_1_str.find_first_of(delimiters, current) != string::npos) {
    extract_res_nums(argv_1_str, selected_residues);
    n = selected_residues.size();
  } else {
    n = atoi(argv[1]);
  }

  ifstream fragments_stream(argv[2]);
  ifstream reference_stream(argv[3]);
  string current_line;
  bool binary_mode = false;

  float rmsd_limit = atof(argv[4]);
  bool filter = (rmsd_limit != 0.0);

  if (argc == 6) {
    binary_mode = string(argv[5]) == "-bin";
  }
  if (not fragments_stream.is_open()) {
    cerr << __FILE__ << ": " << __LINE__
         << ": ERROR: can't read file: " << string(argv[2]) << endl;
    exit(1);
  }
  if (not reference_stream.is_open()) {
    cerr << __FILE__ << ": " << __LINE__
         << ": ERROR: can't read file: " << string(argv[3]) << endl;
    exit(1);
  }

  // load the reference fragment
  while (getline(reference_stream, current_line)) {
    ref_frag.push_back(current_line);
    ++m;
  }
  reference_stream.close();
  if (selected_residues.size() > 0) {
    vector<string> new_ref_frag;
    // cout << "ref BEFORE reduce" << endl;
    // for (size_t i = 0 ; i < ref_frag.size() ; ++i) {
    //   cout << ref_frag[i] << endl;
    // }
    reduce(ref_frag, new_ref_frag, selected_residues);
    ref_frag = new_ref_frag;
    // cout << "ref AFTER reduce" << endl;
    // for (size_t i = 0 ; i < ref_frag.size() ; ++i) {
    //   cout << ref_frag[i] << endl;
    // }
  }
  // - 1 for header line with frag id, / 4 because 4 bb atoms per fragment
  int nb_read = (ref_frag.size() - 1) / 4;
  if (nb_read != n) {
    cerr << __FILE__ << ": " << __LINE__
         << ": FATAL: reference fragment has " << nb_read
         << " residues but expected: " << n << endl;
    exit(1);
  }

  float rms;
  _mLen  = 4 * n;
  _mPOS1 = MatInit(3, _mLen);
  _mPOS2 = MatInit(3, _mLen);
  SimPDB reference(ref_frag, empty_coords, _mLen, _mPOS1);
  int i = 0;

  if (binary_mode) {

    if (selected_residues.size() > 0) {
      cerr << __FILE__ << ": " << __LINE__
           << ": FATAL: holes in fragments in binary mode not yet \
                   supported" << endl;
      exit(1);
    }

    uint32_t raw_float;
    float x;
    string header;
    vector<unsigned char> payload;
    int read_size = m * 3 * 4; // #atoms * #dimensions * #bytes per coord.
    bool can_read = false;

    do {
      can_read = raw_read_one_binary_fragment(read_size,
                                              fragments_stream,
                                              header,
                                              payload);
      //cerr << "read: " << header << endl;
      vector<string> curr_frag;
      curr_frag.push_back(header);
      vector<double> curr_coords;
      curr_coords.reserve(m * 3);

      for (int j = 0 ; j < read_size ; ) { // decode binary coordinates

        raw_float = payload[j++];
        raw_float = (raw_float << 8) + payload[j++];
        raw_float = (raw_float << 8) + payload[j++];
        raw_float = (raw_float << 8) + payload[j++];

        memcpy(&x, &raw_float, 4);
        //cerr << "x: " << x << endl;
        curr_coords.push_back(x);
      }

      SimPDB sim(curr_frag, curr_coords, _mLen, _mPOS2);
      // avoid very small floats (almost 0.0) being printed out in scientific
      // notation while the rest is not
      if (filter) {
        // if (lower_bound_rmsd() <= rmsd_limit) {
        rms = rmsd(0, i);
        if (rms <= rmsd_limit) {
          printf("%.6f %s\n",
                 rms,
                 // skip "TER " at beginning of fragment identifier line
                 curr_frag[0].substr(4).c_str());
        }
        // }
      } else {
        printf("%.6f %s\n",
               rmsd(0, i),
               // skip "TER " at beginning of fragment identifier line
               curr_frag[0].substr(4).c_str());
      }
      ++i;

    } while (can_read);

  } else {
    while (getline(fragments_stream, current_line)) { // read fragment id

      vector<string> curr_frag;
      curr_frag.push_back(current_line);

      for (int j = 0 ; j < m ; ++j) { // read fragment atoms
        getline(fragments_stream, current_line);
        curr_frag.push_back(current_line);
      }
      if (selected_residues.size() > 0) {
        vector<string> new_curr_frag;
        // cout << "curr BEFORE reduce" << endl;
        // for (size_t i = 0 ; i < curr_frag.size() ; ++i) {
        //   cout << curr_frag[i] << endl;
        // }
        reduce(curr_frag, new_curr_frag, selected_residues);
        curr_frag = new_curr_frag;
        // cout << "curr AFTER reduce" << endl;
        // for (size_t i = 0 ; i < curr_frag.size() ; ++i) {
        //   cout << curr_frag[i] << endl;
        // }
      }

      SimPDB sim(curr_frag, empty_coords, _mLen, _mPOS2);

      if (filter) {
        // if (lower_bound_rmsd() <= rmsd_limit) {
        rms = rmsd(0, i);
        if (rms <= rmsd_limit) {
          printf("%.6f %s\n",
                 rms,
                 // skip "TER " at beginning of fragment identifier line
                 curr_frag[0].substr(4).c_str());
        }
        // }
      } else {
        printf("%.6f %s\n",
               rmsd(0, i),
               // skip "TER " at beginning of fragment identifier line
               curr_frag[0].substr(4).c_str());
      }
      ++i;
    }
  }
  fragments_stream.close();

  MatDestroy(&_mPOS1);
  MatDestroy(&_mPOS2);

  return 0;
}
