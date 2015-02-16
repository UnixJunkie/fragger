#!/bin/bash

# only keep backbone atoms
# usage: bb_filter.sh out_dir in.pdb

OUT=$1/`basename $2`

TMP=`mktemp`
trap "rm -f $TMP" EXIT

egrep '^ATOM  ' $2 > $TMP

pdbset xyzin $TMP xyzout $OUT 1>/dev/null << EOF
EXCLUDE SIDE
EXCLUDE WATer
EXCLUDE HYDROGENS
EXCLUDE HEADERS
PICK N CA C O
RENUMBER 1
EOF
