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

module L  = List
module MU = My_utils
module S  = String

(* to get a pdb_chain from disk, kind of inode *)
type pdb_chain_desc =
    { pdb_chain      : string ; (* identifier *)
      offset         : int64  ; (* its start in the DB file *)
      ter_line_len   : int    ;
      atom_lines_len : int    ;
      sequence       : string ; (* amino acid one letter codes *)
      bin_coords     : string }

(* a pdb_chain once read from disk *)
type pdb_chain = string array

let block_length bdesc =
  bdesc.ter_line_len + bdesc.atom_lines_len

let compare_frag_offset f1 f2 =
  Int64.ascending f1.offset f2.offset

let deconstruct_id id = match S.split ~on:'_' id with
  | [pdb; chain] -> (pdb, chain.[0], "c")
  | [pdb; chain; i] -> (pdb, chain.[0], i)
  | _ -> failwith ("indexed_frag.ml: deconstruct_id: invalid id: " ^ id)

(* "TER 1L2Y_A ABCD..." -> ("1L2Y", 'A', "ABCD...") *)
let deconstruct_ter_line l = match S.split l ~on:' ' with
  | [_ter; pdb_chain; seq] ->
    let pdb, chain, maybe_index = deconstruct_id pdb_chain in
    if maybe_index = "c" then
      (pdb, chain, seq)
    else
      failwith ("indexed_frag.ml: deconstruct_ter_line: \
                 we should not have an index at this stage: " ^ l)
  | _ -> failwith ("indexed_frag.ml: deconstruct_ter_line: invalid line: " ^ l)

(* which fragments to generate *)
type frag_id = Chain | All | Any | Just of int

let lookup index id =
  let pdb, chain, maybe_index = deconstruct_id id in
  let m_index = match maybe_index with
    | "c" -> Chain
    | "*" -> All
    | "~" -> Any
    | i ->
      try
        if 4 <> S.length i then
          raise (Failure ""); (* fragment indexes have 4 digits *)
        let requested = Int.of_string i in
        if requested > 0
        then Just requested
        else raise (Failure "")
      with Failure _ ->
        failwith ("indexed_frag.ml: invalid fragment index: " ^ i) in
  let pdb_chain = sprintf "%s_%c" pdb chain in
  let desc = MU.find Fn.id index pdb_chain in
  (desc, m_index)

(* any fragment sequence is OK *)
let no_sequence_filter _ = true

(* manage the string passed by users via the -seq option *)
let sequence_filter filt = match filt with
  | "" -> no_sequence_filter
  | _  ->
    let reg = Str.regexp filt in
    (fun seq -> Str.string_match reg seq 0)
