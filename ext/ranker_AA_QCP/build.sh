#!/bin/bash

set -x

g++ -o SimpPDB.o -c -g -O3 -W -Wall src/SimpPDB.cc
g++ -o qcprot.o  -c -g -O3 -W -Wall src/qcprot.cc
g++ -o ranker.o  -c -g -O3 -W -Wall src/ranker.cc
g++ -o ranker_aa ranker.o SimpPDB.o qcprot.o

# # profile
# g++ -o SimpPDB.o -pg -c -g -O3 -W -Wall src/SimpPDB.cc
# g++ -o qcprot.o  -pg -c -g -O3 -W -Wall src/qcprot.cc
# g++ -o ranker.o  -pg -c -g -O3 -W -Wall src/ranker.cc
# g++ -o ranker_aa -pg ranker.o SimpPDB.o qcprot.o
