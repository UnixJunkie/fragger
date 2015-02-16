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
open Printf
open Indexed_frag
open Log

module A      = Array
module AA     = Amino_acid
module Buff   = Buffer
module F      = Filename
module HT     = Caml.Hashtbl
module L      = List
module MU     = My_utils
module Opt    = Option
module PDB    = Pdb_parser
module S      = String
module StrSet = String.Set

(* how we read and write atom coordinates *)
type io_mode = Txt2txt | Bin2txt | Bin2bin

let io_mode_of_string = function
  | "t2t" -> Txt2txt
  | "b2t" -> Bin2txt
  | "b2b" -> Bin2bin
  | m -> failwith ("frag_db.ml: unsupported io_mode: " ^ m)

(* x, y, z: 3 floats * 4 bytes per float = 12 bytes *)
let atom_len = 12;;

(* counter of output fragments *)
let count = ref 0;;

let get_count () =
  !count

let reset_count () =
  count := 0

(* retrieve a pdb_chain from the database
   !!! PERFORMANCE CRITICAL FUNCTION !!! *)
let get_pdb_chain frag_db_fd frag =
  let block_len = block_length frag in
  let frag_str = String.create block_len in
  (* seek it *)
  ignore(Unix.lseek frag_db_fd frag.offset Unix.SEEK_SET);
  (* read it *)
  let was_read = MU.really_read frag_db_fd frag_str block_len in
  if was_read <> block_len then begin
    fatal (lazy (
      sprintf "frag_db.ml: get_pdb_chain: could not read whole fragment: %s"
        frag.pdb_chain));
    exit 1;
  end;
  let atom_line_len = 81 in
  let n_atom_lines = frag.atom_lines_len / atom_line_len in
  let atom_lines = A.create n_atom_lines "" in
  (* put ATOM lines into an array *)
  for i = 0 to n_atom_lines - 1 do
    atom_lines.(i) <-
      S.sub frag_str (frag.ter_line_len + (atom_line_len * i)) atom_line_len;
  done;
  atom_lines

let fragment_sequence frag nb_res index =
  S.sub frag.sequence (index - 1) nb_res

(* fragments indexes start at 1 (as the first residue number in a PDB) *)
let output_one_frag_priv no_coords frag seq_filter out pdb_chain nb_res index =
  if nb_res = 0
  then failwith "frag_db.ml: output_one_frag_priv: 0 length fragment";
  let sub_seq = fragment_sequence frag nb_res index in
  if seq_filter sub_seq then begin
    incr count;
    (* TER line *)
    if no_coords then begin
      fprintf out "TER %s_%04d %s\n" frag.pdb_chain index sub_seq
    end else begin
      fprintf out "TER %s_%04d\n" frag.pdb_chain index;
      (* ATOM lines *)
      (*           4 backbone atoms per residue *)
      let start  = 4 * (index - 1)          in
      let finish = start + (4 * nb_res) - 1 in
      for i = start to finish do
        fprintf out "%s" pdb_chain.(i)
      done
    end
  end

(* possible index range for some given fragment lines and a fragment length *)
let max_index frag nb_res =
  (S.length frag.sequence) - (nb_res - 1)

let create_pdb_line atom_num atom_name res_name chain res_num x y z =
(* example:
"ATOM   2372  CA  PHE B 307     117.704  -5.451 100.082  1.00  6.56\
      B    C  " *)
  sprintf
    "ATOM  %5d%s  %s %c%4d    %8.3f%8.3f%8.3f  1.00 30.00              \n"
    atom_num atom_name res_name chain res_num x y z

let output_a_frag no_coords pdb_chain frag seq_filter out nb_res index =
  let max_i = max_index frag nb_res in
  if index < 1 || index > max_i then
    warn (lazy (
      sprintf "frag_db.ml: output_a_frag of len %d: impossible: %s %d"
        nb_res frag.pdb_chain index))
  else
    output_one_frag_priv no_coords frag seq_filter out pdb_chain nb_res index

let output_a_frag_b2t no_coords frag seq_filter out nb_res index =
  let sub_seq = fragment_sequence frag nb_res index in
  if seq_filter sub_seq then begin
    incr count;
    let _pdb, chain, _index = deconstruct_id frag.pdb_chain in
    (* TER line *)
    if no_coords then
      fprintf out "TER %s_%04d %s\n" frag.pdb_chain index sub_seq
    else (
      (* for speed, we ignore the fragment sequence in this output format *)
      fprintf out "TER %s_%04d\n" frag.pdb_chain index;
      (* ATOM lines *)
      (*          4 backbone atoms per residue *)
      let start = 4 * (index - 1) in
      let finish = start + (4 * nb_res) - 1 in
      for i = start to finish do
        let coords = S.sub frag.bin_coords (atom_len * i) atom_len in
        let x_bin = S.sub coords 0 4 in
        let y_bin = S.sub coords 4 4 in
        let z_bin = S.sub coords 8 4 in
        let x = MU.float32_of_str x_bin in
        let y = MU.float32_of_str y_bin in
        let z = MU.float32_of_str z_bin in
        let atom_name = match i mod 4 with
          | 0 -> "  N "
          | 1 -> "  CA"
          | 2 -> "  C "
          | 3 -> "  O "
          | _ -> failwith "output_a_frag_b2t: IMPOSSIBLE"
        in
        let res_num = (i / 4) + 1 in
        let res = frag.sequence.[res_num - 1] in
        let res_name =
          AA.string_of_aa_three
            (AA.aa_three_of_string
               (String.make 1 res)) in
        let atom_line =
          create_pdb_line (i + 1) atom_name res_name chain res_num x y z in
        fprintf out "%s" atom_line
      done)
    end

let output_a_frag_b2b no_coords frag seq_filter out nb_res index =
  let sub_seq = fragment_sequence frag nb_res index in
  if seq_filter sub_seq then begin
    incr count;
    (* TER line *)
    if no_coords then
      fprintf out "TER %s_%04d %s\n" frag.pdb_chain index sub_seq
    else
      (* for speed, we ignore the fragment sequence in this output format *)
      fprintf out "TER %s_%04d\n" frag.pdb_chain index;
      (* ATOM lines binary coordinates *)
      (*          4 backbone atoms per residue *)
      let start = 4 * atom_len * (index - 1) in (* indexes start at 1 *)
      let len = 4 * atom_len * nb_res in
      let coords = S.sub frag.bin_coords start len in
      fprintf out "%s" coords
  end

let output_any_frag no_coords frag seq_filter out pdb_chain nb_res =
  let max_i = max_index frag nb_res in
  if max_i < 1 then
    warn (lazy (sprintf "skipped %s (too short)" frag.pdb_chain))
  else
    (* any is in [1 ; max_index] *)
    let any = 1 + Random.int max_i in
    output_one_frag_priv no_coords frag seq_filter out pdb_chain nb_res any

let output_chain no_coords frag out pdb_chain =
  incr count;
  (* TER line *)
  if no_coords then
    fprintf out "TER %s %s\n" frag.pdb_chain frag.sequence
  else begin
    fprintf out "TER %s\n" frag.pdb_chain;
    (* ATOM lines *)
    A.iter ~f:(fprintf out "%s") pdb_chain
  end

let output_all_frags_t2t no_coords frag seq_filter out pdb_chain nb_res =
  let max_i = max_index frag nb_res in
  for i = 1 to max_i do
    output_one_frag_priv no_coords frag seq_filter out pdb_chain nb_res i
  done

let output_all_frags_b2t no_coords frag seq_filter out nb_res =
  let max_i = max_index frag nb_res in
  for i = 1 to max_i do
    output_a_frag_b2t no_coords frag seq_filter out nb_res i
  done

let output_all_frags_b2b no_coords frag seq_filter out nb_res =
  let max_i = max_index frag nb_res in
  for i = 1 to max_i do
    output_a_frag_b2b no_coords frag seq_filter out nb_res i
  done

let output_frags
    no_coords mode seq_filter out frag_db_fd nb_res frag maybe_ids =
  let mode_specialized, m_pdb_chain = match mode with
    | Txt2txt -> let pdb_chain = get_pdb_chain frag_db_fd frag in
                 (output_a_frag no_coords pdb_chain frag seq_filter out nb_res,
                  Some pdb_chain)
    | Bin2txt -> output_a_frag_b2t no_coords frag seq_filter out nb_res, None
    | Bin2bin -> output_a_frag_b2b no_coords frag seq_filter out nb_res, None
  in
  L.iter maybe_ids
    ~f:(function
        | Just i -> mode_specialized i
        | Any -> let pdb_chain = Opt.value_exn m_pdb_chain in
                 output_any_frag no_coords frag seq_filter out pdb_chain nb_res
        | All -> (match mode with
            | Txt2txt ->
              let pdb_chain = Opt.value_exn m_pdb_chain in
              output_all_frags_t2t
                no_coords frag seq_filter out pdb_chain nb_res
            | Bin2txt ->
              output_all_frags_b2t no_coords frag seq_filter out nb_res
            | Bin2bin ->
              output_all_frags_b2b no_coords frag seq_filter out nb_res)
        | Chain -> let pdb_chain = Opt.value_exn m_pdb_chain in
                   output_chain no_coords frag out pdb_chain)

let output_several_frags
    no_coords mode seq_filter index frag_db_fd nb_res requested out =
  (* varies each run *)
  let _ = Random.self_init() in
  (* lookup all requested and store them in an array *)
  let table = A.of_list (L.rev_map requested ~f:(lookup index)) in
  (* sort by offset *)
  A.sort table
    ~cmp:(fun (pdb_chain1, _) (pdb_chain2, _) ->
            compare_frag_offset pdb_chain1 pdb_chain2);
  (* group by same offset/pdb_chain *)
  let last =
    A.fold table ~init:None
      ~f:(fun acc (pdb_chain, id) -> match acc with
            | None -> Some (pdb_chain, [id])
            | Some (prev_pdb_chain, ids) ->
              if Int64.equal pdb_chain.offset prev_pdb_chain.offset then
                Some (prev_pdb_chain, id :: ids)
              else
                let _ =
                  output_frags no_coords mode seq_filter out frag_db_fd nb_res
                    prev_pdb_chain ids in
                Some (pdb_chain, [id]))
  in
  match last with
    | None -> ()
    | Some (prev_pdb_chain, ids) ->
      output_frags
        no_coords mode seq_filter out frag_db_fd nb_res prev_pdb_chain ids

let output_one_frag
    no_coords mode seq_filter index frag_db_fd nb_res requested out =
  let pdb_chain, m_index = lookup index requested in
  output_frags
    no_coords mode seq_filter out frag_db_fd nb_res pdb_chain [m_index]

(*
  let frag_regexp = Str.regexp "^TER " (* a fragment identifier line *)
  let is_frag_line l = Str.string_match frag_regexp l 0
*)

let load_db_index db_file debug =
  let debug_f index =
    (* print out content of the fragments index as text *)
    HT.iter
      (fun _id desc ->
        printf "frag:%s start:%Ld len:%d\n"
          desc.pdb_chain
          desc.offset
          (block_length desc))
      index
  in
  let bin_idx_file = db_file ^ ".idx" in
  if Sys.file_exists bin_idx_file = `Yes then begin
    info (lazy "loading DB index from previous run ...");
    let frags_index =
      In_channel.with_file ~binary:true bin_idx_file
        ~f:(fun input -> Marshal.from_channel input) in
    info (lazy "loaded");
    MU.guarded_tap debug debug_f frags_index
  end else begin
    failwith ("Recreate your DB: could not open index: " ^ bin_idx_file);
  end
