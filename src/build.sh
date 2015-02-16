#!/bin/bash

#set -x

rm -f setup.*
oasis setup
ocaml setup.ml -configure
make
