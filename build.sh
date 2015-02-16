#!/bin/bash

set -x

# OCaml libraries we depend on
opam update
opam install core batteries parmap dolog oasis

# our fast RMSD computation tool
cd ext/ranker_AA_QCP/
./build.sh

# Fragger
cd ../../src
rm -rf _build
./build.sh
