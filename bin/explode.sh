#!/bin/bash

if [ "$#" -ne "2" ] ; then
    echo "usage: "$0" fragments_file output_directory"
    exit 1
fi

awk -v d=$2 ' $1=="TER" {file=sprintf("%s/%s.pdb",d,$2); print $0 > file;} $1=="ATOM" {print $0 >> file;}' $1
