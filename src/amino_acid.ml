(* Copyright (C) 2013, Zhang Initiative Research Unit,
 * Advance Science Institute, Riken
 * 2-1 Hirosawa, Wako, Saitama 351-0198, Japan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License, with
 * the special exception on linking described in file LICENSE.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. *)

open Core.Std

type aa_three =
  | ALA
  | ARG
  | ASN
  | ASP
  | CYS
  | GLU
  | GLN
  | GLY
  | HIS
  | ILE
  | LEU
  | LYS
  | MET
  | PHE
  | PRO
  | SER
  | THR
  | TRP
  | TYR
  | VAL

type aa_one =
  | A
  | R
  | N
  | D
  | C
  | E
  | Q
  | G
  | H
  | I
  | L
  | K
  | M
  | F
  | P
  | S
  | T
  | W
  | Y
  | V

let three_to_one = function
  | ALA -> A
  | ARG -> R
  | ASN -> N
  | ASP -> D
  | CYS -> C
  | GLU -> E
  | GLN -> Q
  | GLY -> G
  | HIS -> H
  | ILE -> I
  | LEU -> L
  | LYS -> K
  | MET -> M
  | PHE -> F
  | PRO -> P
  | SER -> S
  | THR -> T
  | TRP -> W
  | TYR -> Y
  | VAL -> V

let one_to_three = function
  | A -> ALA
  | R -> ARG
  | N -> ASN
  | D -> ASP
  | C -> CYS
  | E -> GLU
  | Q -> GLN
  | G -> GLY
  | H -> HIS
  | I -> ILE
  | L -> LEU
  | K -> LYS
  | M -> MET
  | F -> PHE
  | P -> PRO
  | S -> SER
  | T -> THR
  | W -> TRP
  | Y -> TYR
  | V -> VAL

exception Unknown_aa_three of string

let aa_three_of_string = function
  | "ALA" -> ALA
  | "ARG" -> ARG
  | "ASN" -> ASN
  | "ASP" -> ASP
  | "CYS" -> CYS
  | "GLU" -> GLU
  | "GLN" -> GLN
  | "GLY" -> GLY
  | "HIS" -> HIS
  | "ILE" -> ILE
  | "LEU" -> LEU
  | "LYS" -> LYS
  | "MET" -> MET
  | "PHE" -> PHE
  | "PRO" -> PRO
  | "SER" -> SER
  | "THR" -> THR
  | "TRP" -> TRP
  | "TYR" -> TYR
  | "VAL" -> VAL
  | "A" -> ALA
  | "R" -> ARG
  | "N" -> ASN
  | "D" -> ASP
  | "C" -> CYS
  | "E" -> GLU
  | "Q" -> GLN
  | "G" -> GLY
  | "H" -> HIS
  | "I" -> ILE
  | "L" -> LEU
  | "K" -> LYS
  | "M" -> MET
  | "F" -> PHE
  | "P" -> PRO
  | "S" -> SER
  | "T" -> THR
  | "W" -> TRP
  | "Y" -> TYR
  | "V" -> VAL
  | unk   -> raise (Unknown_aa_three unk)

exception Unknown_aa_one of string

let aa_one_of_string = function
  | "ALA" -> A
  | "ARG" -> R
  | "ASN" -> N
  | "ASP" -> D
  | "CYS" -> C
  | "GLU" -> E
  | "GLN" -> Q
  | "GLY" -> G
  | "HIS" -> H
  | "ILE" -> I
  | "LEU" -> L
  | "LYS" -> K
  | "MET" -> M
  | "PHE" -> F
  | "PRO" -> P
  | "SER" -> S
  | "THR" -> T
  | "TRP" -> W
  | "TYR" -> Y
  | "VAL" -> V
  | "A" -> A
  | "R" -> R
  | "N" -> N
  | "D" -> D
  | "C" -> C
  | "E" -> E
  | "Q" -> Q
  | "G" -> G
  | "H" -> H
  | "I" -> I
  | "L" -> L
  | "K" -> K
  | "M" -> M
  | "F" -> F
  | "P" -> P
  | "S" -> S
  | "T" -> T
  | "W" -> W
  | "Y" -> Y
  | "V" -> V
  | unk -> raise (Unknown_aa_one unk)

let string_of_aa_three = function
  | ALA -> "ALA"
  | ARG -> "ARG"
  | ASN -> "ASN"
  | ASP -> "ASP"
  | CYS -> "CYS"
  | GLU -> "GLU"
  | GLN -> "GLN"
  | GLY -> "GLY"
  | HIS -> "HIS"
  | ILE -> "ILE"
  | LEU -> "LEU"
  | LYS -> "LYS"
  | MET -> "MET"
  | PHE -> "PHE"
  | PRO -> "PRO"
  | SER -> "SER"
  | THR -> "THR"
  | TRP -> "TRP"
  | TYR -> "TYR"
  | VAL -> "VAL"

let string_of_aa_one = function
  | A -> "A"
  | R -> "R"
  | N -> "N"
  | D -> "D"
  | C -> "C"
  | E -> "E"
  | Q -> "Q"
  | G -> "G"
  | H -> "H"
  | I -> "I"
  | L -> "L"
  | K -> "K"
  | M -> "M"
  | F -> "F"
  | P -> "P"
  | S -> "S"
  | T -> "T"
  | W -> "W"
  | Y -> "Y"
  | V -> "V"

let char_of_aa_one = function
  | A -> 'A'
  | R -> 'R'
  | N -> 'N'
  | D -> 'D'
  | C -> 'C'
  | E -> 'E'
  | Q -> 'Q'
  | G -> 'G'
  | H -> 'H'
  | I -> 'I'
  | L -> 'L'
  | K -> 'K'
  | M -> 'M'
  | F -> 'F'
  | P -> 'P'
  | S -> 'S'
  | T -> 'T'
  | W -> 'W'
  | Y -> 'Y'
  | V -> 'V'
