#!/bin/bash

if [ "$#" -ne "4" ]
then
    echo "usage: three2one.sh fragments_file fragment_length window_number \
output_file"
    exit 1
fi

awk -v m=$2 -v r=$3 '
function three2one(residue)
{
    switch(residue)
    {
        case "ALA": return "A";break;
        case "ARG": return "R";break;
        case "ASN": return "N";break;
        case "ASP": return "D";break;
        case "CYS": return "C";break;
        case "GLU": return "E";break;
        case "GLN": return "Q";break;
        case "GLY": return "G";break;
        case "HIS": return "H";break;
        case "ILE": return "I";break;
        case "LEU": return "L";break;
        case "LYS": return "K";break;
        case "MET": return "M";break;
        case "PHE": return "F";break;
        case "PRO": return "P";break;
        case "SER": return "S";break;
        case "THR": return "T";break;
        case "TRP": return "W";break;
        case "TYR": return "Y";break;
        case "VAL": return "V";break;
        default: print "No matching residue:", residue; exit 1;
    }
}

BEGIN{ind=1;}
$1=="REMARK" {rmsd=$2; printf("%d %f ", r, $2);}
$1=="ATOM" && $3=="N" { residue = substr($0,18,3); 
if (ind<m) { ind+=1; printf("%c ", three2one(residue));} 
else {ind=1;printf("%c\n", three2one(residue));}} ' $1 > $4
